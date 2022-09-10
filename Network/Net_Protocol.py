from .Packet import Net_Packet
from .Address import *
from Sim.Simulator import *

class Net_Protocol:
    
    def __init__(self):
        self.node = None
        self.routing = None

    def set_node(self, node):
        self.node = node
        node.net_protocol = self

    def set_routing(self, routing):
        self.routing = routing
        routing.address = self.node.device.address

    def send_data(self, source, destination, data_packet):
        np_header = Net_Packet.Header(source, destination)
        if destination == Addresses.BROADCAST:
            net_packet = Net_Packet(data_packet)
            self.node.net_device.send(net_packet, destination, None)
            return
        net_route = self.routing.route_output(np_header)
        if net_route:
            self.send_out(net_route, np_header, data_packet)
        else:
            print('No route to destination node. Drop.')

    def receive_net_packet(self, net_packet, fr, to):
        self.routing.route_input(
            net_packet, fr, to, self.send_data, self.receive_data, self.forward)

    def receive_data(self, net_packet):
        node = net_packet.header.source.get_node()
        self.node.receive_from(node, net_packet.payload)

    def send_out(self, net_route, np_header, data_packet):
        net_packet = Net_Packet(np_header, data_packet)
        if net_route.next_hop == Addresses.LOOPBACK:
            fr = Addresses.LOOPBACK
            to = Addresses.LOOPBACK
            Simulator.schedule_now(self.receive_data, fr, to, net_packet)
            return True
        else:
            return self.node.net_device.send(net_packet, net_route.next_hop, net_route.tx_params)

    def forward(self, net_route, net_packet):
        np_header = net_packet.header
        data_packet = net_packet.payload
        return self.send_out(net_route, np_header, data_packet)

    class Helper:

        def __init__(self):
            self.routing_helper = None

        def set_routing_helper(self, routing_helper):
            self.routing_helper = routing_helper

        def install(self, nodes):
            for node in nodes.get_nodes:
                net_protocol = Net_Protocol()
                net_protocol.routing = self.routing_helper.create()
                net_protocol.set_node(node)

class Net_Route:
    def __init__(self):
        self.source = None
        self.destination = None
        self.next_hop = None
        self.tx_params = None
