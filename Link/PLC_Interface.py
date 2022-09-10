from Sim.Simpy_Core import *
from Sim.Simulator import *
from Network.Packet import *
from Sim.Time import *
from Sim.Logging import *
import numpy as np


class PLC_Interface:
    def __init__(self):
        self.node = None
        self.mac = None
        self.phy = None
        self.channel = None
        self.state = self.State.IDLE
        self.tx_power = 20
        self.tx_psd = None
        self.tx_signal = None
        # cca threshold, -62dBm
        self.cca_threshold = 6.3e-7
        self.cca_neighbours = []
        self.rx_paras = []
        self.rx_log = []
        PLC_Interfaces.add(self)

    def set_node(self, node):
        self.node = node

    def set_mac(self, mac):
        self.mac = mac

    def set_phy(self, phy):
        self.phy = phy

    def set_channel(self, channel):
        self.channel = channel

    def set_tx_power(self, tx_power):
        self.tx_power = tx_power

    def update_cca_neighbors(self):
        channel = self.channel
        tx_psd = self.phy.robo_psd_alloc(self.tx_power)
        self.cca_neighbours = channel.transmission_range_calc(
            self, tx_psd, self.cca_threshold)

    def cca(self):
        cca_result = 1
        for interface in self.cca_neighbours:
            if interface.state != self.State.IDLE:
                cca_result = 0
                break
        return cca_result

    def send(self, mac_frame, tx_para, callback=None):
        phy_frame = Phy_Frame(mac_frame)
        self.tx_init(phy_frame, tx_para)
        self.tx_start(phy_frame, callback)

    def rx_para_calc(self, tx_psd, rx_psd):
        signals = PLC_Signal.signals
        channel = self.channel
        interfer_psd = np.zeros(len(rx_psd))
        for signal in signals:
            if self not in signal.interface.cca_neighbours:
                interfer_rx_psd = channel.transmission_calc(
                    signal.interface, self, signal.tx_psd)
                interfer_psd += interfer_rx_psd
        est_para, rx_sinr = self.phy.channel_est(
            tx_psd, rx_psd, interfer_psd)
        self.rx_paras.append(self.Rx_Para(est_para, rx_sinr))

    def tx_init(self, phy_frame, tx_para):
        if tx_para.tx_mode == self.Tx_Mode.Broadcast:
            self.tx_psd = self.phy.broadcast_mode(
                phy_frame.preamble, self.tx_power)
        else:
            self.tx_psd = self.phy.unitcast_mode(
                phy_frame.preamble, self.tx_power, tx_para.to_para_est)

    def tx_start(self, phy_frame, callback=None):
        self.state = self.State.TX
        self.tx_signal = PLC_Signal(self, self.tx_psd)
        channel = self.channel
        PLC_Signal.add(self.tx_signal)
        tx_duration, num_blocks = self.phy.tx_duration_calc(phy_frame)
        for interface in self.cca_neighbours:
            rx_psd = channel.transmission_calc(
                self, interface, self.tx_psd)
            interface.rx_start(self.tx_psd, rx_psd, tx_duration, num_blocks)
        Simulator.schedule(tx_duration.sim_time,
                           self.tx_end, phy_frame, callback)

    def tx_end(self, phy_frame, callback=None):
        PLC_Signal.delete(self.tx_signal)
        self.state = self.State.IDLE
        for interface in self.cca_neighbours:
            interface.rx_end(phy_frame)
        if callback:
            callback.succeed()

    def rx_start(self, tx_psd, rx_psd, rx_duration, num_blocks):
        self.state = self.State.RX
        self.rx_paras = []
        for i in range(num_blocks):
            slot = rx_duration.sim_time/num_blocks
            Simulator.schedule(i*slot, self.rx_para_calc, tx_psd, rx_psd)

    def rx_end(self, phy_frame):
        self.state = self.State.IDLE
        rx_sinrs = [para.rx_sinr for para in self.rx_paras]
        prob_rx_success = self.phy.error_rate_model.packet_success_rate_calc(
            phy_frame, rx_sinrs)
        rx_log = Logging.Rx_Log(self, prob_rx_success,
                                copy.deepcopy(rx_sinrs[0]))
        self.rx_log.append(rx_log)
        print('prob rx', prob_rx_success)
        if random.random() <= prob_rx_success:
            print('node', self.node.id, 'plc interface received at ',
                  Simulator.now().value, Simulator.now().type)

            if self.mac:
                rx_para = copy.deepcopy(self.rx_paras)
                self.receive(phy_frame, rx_para)
        else:
            print('node', self.node.id, 'plc interface decode failed at ',
                  Simulator.now().value, Simulator.now().type)

    def receive(self, phy_frame):
        self.mac.receive(phy_frame.payload)

    class Tx_Mode:
        Broadcast = 0
        Unicast = 1

    class Tx_Para:
        def __init__(self, tx_mode, to_para_est=None):
            self.tx_mode = tx_mode
            self.to_para_est = to_para_est

    class Rx_Para:
        def __init__(self, est_para, rx_sinr):
            self.est_para = est_para
            self.rx_sinr = rx_sinr

    class State:
        TX = 0
        RX = 1
        IDLE = 2


class PLC_Interfaces:
    interfaces = []

    @classmethod
    def add(cls, interface):
        cls.interfaces.append(interface)

    @classmethod
    def get(cls, node):
        for interface in cls.interfaces:
            if interface.node == node:
                return interface
        return None

    @classmethod
    def get_all(cls):
        return cls.interfaces


class PLC_Signal:
    signals = []

    def __init__(self, interface, tx_psd):
        self.interface = interface
        self.tx_psd = tx_psd

    @classmethod
    def get_all(cls):
        return cls.signals

    @classmethod
    def add(cls, signal):
        if signal not in cls.signals:
            cls.signals.append(signal)

    @classmethod
    def delete(cls, signal):
        if signal in cls.signals:
            cls.signals.remove(signal)
