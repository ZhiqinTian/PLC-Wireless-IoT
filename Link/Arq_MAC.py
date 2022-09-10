from Sim.Simpy_Core import *
from Sim.Simulator import *
from Network.Packet import *
from Network.Address import Addresses
from Sim.Time import *
import copy


class Arq_MAC:
    def __init__(self):
        self.net_device = None
        self.phy_interface = None
        self.neighbors = {}
        self.tx_queue = Queue()
        self.rx_queue = Queue()
        self.max_rx_caches = 10
        self.max_retries = 5
        self.max_timeout = MilliSecond(1)
        self.slot = MicroSecond(10)
        self.ack_event = None
        self.tx_callback = None
        self.ack_sending = None

    def set_net_device(self, device):
        self.net_device = device

    def set_interface(self, interface):
        self.phy_interface = interface
        interface.mac = self

    def send(self, destination, net_packet):
        if self.net_device:
            source = self.net_device.address
        else:
            source = None
        mac_header = MAC_Frame.Header(source, destination)
        mac_frame = MAC_Frame(mac_header, net_packet)
        self.tx_queue.append(mac_frame)
        if len(self.tx_queue) == 1:
            Simulator.schedule_now(self.queue_process)

    def sendb_ack(self, mac_frame):
        self.ack_sending = Simulator.create_event()
        k = 1
        while True:
            if self.phy_interface.cca():
                break
            wait_time = MilliSecond(
                self.slot.sim_time*random.randrange(1, k+1))
            yield Simulator.time_delay(wait_time.sim_time)
            k = min(k*2, 16)
        yield Simulator.time_delay(2*self.slot.sim_time)
        tx_mode = self.phy_interface.Tx_Mode.Unicast
        self.phy_interface.send(mac_frame, tx_mode, self.ack_sending)

    def queue_process(self):
        retries = 0
        while self.tx_queue:
            if retries == self.max_retries:
                retries = 0
                dst = self.tx_queue[0].header.destination
                del self.neighbors[dst]
                self.tx_queue.popleft()
            mac_frame = self.tx_queue[0]

            yield Simulator.time_delay(3*self.slot.sim_time)

            # persistent process with exponential backoff
            k = 1
            while True:
                if self.phy_interface.cca():
                    break
                wait_time = MilliSecond(
                    self.slot.sim_time*random.randrange(1, k+1))
                yield Simulator.time_delay(wait_time.sim_time)
                k = min(k*2, 1024)

            if self.ack_sending:
                yield self.ack_sending
                self.ack_sending = None

            self.tx_callback = Simulator.create_event()

            # wait for ack if this is a unicast frame
            dst = mac_frame.header.destination
            if dst != Addresses.BROADCAST:

                tx_mode = self.phy_interface.Tx_Mode.Unicast
                if dst not in self.neighbors:
                    self.neighbors[dst] = self.Link_Est()
                rx_para_est = copy.deepcopy(self.neighbors[dst].to_para)
                tx_para = self.phy_interface.Tx_Para(tx_mode, rx_para_est)
                self.phy_interface.send(mac_frame, tx_para, self.tx_callback)

                self.ack_event = Simulator.create_event()
                self.ack_event.wait_for = mac_frame

                yield AnyOf(Simulator.env, [Simulator.time_delay(self.max_timeout.sim_time), self.ack_event])
                if self.ack_event.triggered:
                    retries = 0
                    self.tx_queue.popleft()
                else:
                    retries += 1
                    backoff = MilliSecond(
                        self.slot.sim_time*random.randrange(1, 2**retries))
                    yield Simulator.time_delay(backoff.sim_time)
            else:
                tx_mode = self.phy_interface.Tx_Mode.Broadcast
                print('node', self.phy_interface.node.id, 'mac frame send at ',
                      Simulator.now().value, Simulator.now().type)
                tx_para = self.phy_interface.Tx_Para(tx_mode)
                self.phy_interface.send(mac_frame, tx_para, self.tx_callback)
                retries = 0
                self.tx_queue.popleft()

            yield self.tx_callback
            self.tx_callback = None
            self.ack_event = None

    def receive_enqueue(self, mac_frame):
        if mac_frame in self.rx_queue:
            return False
        if len(self.rx_queue) == self.max_rx_caches:
            self.rx_queue.popleft()
        self.rx_queue.append(mac_frame)
        return True

    def receive(self, mac_frame, rx_para):
        fr = mac_frame.header.source

        if fr not in self.neighbors:
            self.neighbors[fr] = self.Link_Est()
        self.neighbors[fr].fr_para = rx_para

        destination = mac_frame.header.destination
        if mac_frame.payload.type == 'net_packet':
            print('node', self.phy_interface.node.id, 'mac frame received at ',
                  Simulator.now().value, Simulator.now().type)
            if self.net_device:
                if (destination == Addresses.BROADCAST) or (destination == self.net_device.address):
                    if self.receive_enqueue(self, mac_frame):
                        self.net_device.receive(mac_frame, fr)

                    # ack if this is a unicast frame
                    if destination != Addresses.BROADCAST:
                        mac_header = MAC_Frame.Header(
                            self.net_device.address, fr)
                        mac_para = copy.deepcopy(self.neighbors[fr].fr_para)
                        mac_ack = MAC_Ack(mac_frame, mac_para)
                        mac_frame = MAC_Frame(mac_header, mac_ack)
                        Simulator.schedule_now(self.send_ack, mac_frame)

        elif (mac_frame.payload.type == 'mac_ack') and self.ack_event:
            print('node', self.phy_interface.node.id, 'mac ack received at ',
                  Simulator.now().value, Simulator.now().type)
            if mac_frame.for_frame == self.ack_event.wait_for:
                self.neighbors[fr].to_para = mac_frame.mac_para
                self.ack_event.succeed()

    class Link_Est:
        def __init__(self):
            self.to_para = None
            self.fr_para = None
