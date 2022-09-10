class Net_Device:

    def __init__(self):
        self.node = None
        self.address = None
        self.mac = None
        self.wireless_phy = None
        self.plc_phy = None

    def set_address(self, address):
        self.address = address
        address.net_device = self

    def set_wireless_phy(self, wireless_phy):
        self.wireless_phy = wireless_phy

    def set_plc_phy(self, plc_phy):
        self.plc_phy = plc_phy

    def send(self, net_packet, next_hop, tx_params):
        return self.mac.send(next_hop, net_packet)

    def receive(self, mac_frame):
        net_packet = mac_frame.net_packet
        fr = mac_frame.source
        to = self.address
        self.node.net_protocol.receive_net_packet(net_packet, fr, to)

    def install(self, node):
        self.node = node
        node.net_device = self


class Net_Devices:

    def __init__(self, net_devices=[]):
        self.net_devices = net_devices

    @property
    def get_net_devices(self):
        return self.net_devices

    def get_net_device(self, n):
        return self.net_devices[n]

    def add(self, net_device):
        self.net_devices.append(net_device)

    class Helper:
        def __init__(self):
            self.set_wireless = None
            self.set_plc = None
            self.net_devices = Net_Devices()

        def set_wireless(self, set_wireless):
            self.set_wireless = set_wireless

        def set_plc(self, set_plc):
            self.set_plc = set_plc

        def install(self, nodes):
            for node in nodes:
                net_device = Net_Device()
                set_wireless = self.set_wireless
                set_plc = self.set_plc
                if set_wireless:
                    net_device.wireless_phy = set_wireless.create()
                if set_plc:
                    net_device.plc_phy = set_plc.create()
                net_device.install(node)
                self.net_devices.add(net_device)
