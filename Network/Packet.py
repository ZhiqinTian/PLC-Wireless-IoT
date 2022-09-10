class Byte:
    def __init__(self, value):
        self.type = 'B'
        self.value = value

    @property
    def bits(self):
        return self.value*8

    def to_bit(self):
        return Bit(self.value*8)


class Bit:
    def __init__(self, value):
        self.type = 'b'
        self.value = value

    @property
    def bits(self):
        return self.value

    def to_byte(self):
        return Byte(self.value/8)


class Data_Packet:
    def __init__(self, size=Byte(2000)):
        self.type = 'data_packet'
        self.size = size


class Net_Packet:
    def __init__(self, header=None, payload=None, payload_size=Byte(2000)):
        self.type = 'net_packet'
        if header == None:
            header = self.Header()
        self.header = header
        self.payload = payload
        if payload:
            self.size = Byte(self.header.size.value+self.payload.size.value)
        else:
            self.size = Byte(self.header.size.value+payload_size.value)

    class Header:
        def __init__(self, source=None, destination=None):
            self.source = source
            self.destination = destination
            self.type = 'np_header'
            self.size = Byte(8)


class MAC_Frame:
    def __init__(self, header=None,  payload=None, payload_size=Byte(2000)):
        self.type = 'mac_frame'
        if header == None:
            header = self.Header()
        self.header = header
        self.payload = payload
        if payload:
            self.size = Byte(self.header.size.value+self.payload.size.value)
        else:
            self.size = Byte(self.header.size.value+payload_size.value)

    class Header:
        def __init__(self, source=None, destination=None):
            self.source = source
            self.destination = destination
            self.type = 'mac_header'
            self.size = Byte(9)


class MAC_Ack:
    def __init__(self, for_mac=None, mac_sinr=None):
        self.type = 'mac_ack'
        self.size = Byte(9)
        self.for_mac = for_mac
        self.mac_sinr = mac_sinr


class Phy_Frame:
    def __init__(self, payload=None, payload_size=Byte(2000)):
        self.type = 'phy_frame'
        self.preamble = self.Preamble()
        self.payload = payload

        if payload:
            self.size = self.payload.size
        else:
            self.size = payload_size

    class Preamble:
        def __init__(self):
            self.symbols_num = 5
            self.carriers = None
            self.data_mcs = None
            self.p_mcs = None
