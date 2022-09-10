# define address class
class Address:

    # Constructor: set the broadcast address and the begin of address assignment
    def __init__(self, ad_value):
        self.value = ad_value
        self.net_device = None
        if ad_value!=0xFFFF and ad_value!=0x0000:
            Addresses.add_address(self)

    def __eq__(self, other):
        return self.value == other.value

    # Assignment of address to net-device
    def set_device(self, net_device):
        self.net_device = net_device
        net_device.address = self

    def get_node(self):
        return self.net_device.node

    class Helper:
        def __init__(self, begin):
            self.begin = begin

        def assign(self, net_devices):
            ad_v = self.begin
            for nd in net_devices.get_devices:
                address = Address(ad_v)
                nd.set_address(address)
                ad_v = ad_v+1


# define the broadcast, loopback address and the address collection
class Addresses:
    BROADCAST = Address(0xFFFF)
    LOOPBACK = Address(0x0000)
    addresses = dict()

    @classmethod
    def add_address(cls, address):
        cls.addresses[address.value] = address
    
    @classmethod
    def get_address(cls, address_value):
        return cls.addresses[address_value]
