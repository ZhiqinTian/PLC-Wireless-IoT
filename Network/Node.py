class Node:

    def __init__(self, name=None, ID=None):
        self.name = name
        self.id = ID
        self.net_device = None
        self.net_protocol = None

    def set_net_device(self, net_device):
        self.net_device = net_device
        net_device.node = self

    def set__net_protocol(self, net_protocol):
        self.net_protocol = net_protocol
        net_protocol.node = self

    def send_to(self, node, data_packet):
        source = self.net_device.address
        destination = node.net_device.address
        self.net_protocol.send_data(source, destination, data_packet)

    def receive_from(self,node,data_packet):
        pass

    @classmethod
    def creat_nodes(cls, n):
        nodes = []
        for i in range(n):
            node = Node('n'+str(i+1), i)
            nodes.append(node)
        return nodes
    
    def get_node(cls,nodes,ID):
        for node in nodes:
            if node.id == ID:
                return node
        return None
