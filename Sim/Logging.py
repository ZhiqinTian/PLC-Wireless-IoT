import copy

class Logging:
    phy_log = False
    mac_log = False
    net_log = False
    app_log = False

    @classmethod
    def log_phy(cls, msg):
        if cls.phy_log == True:
            pass
    
    @classmethod
    def log_mac(cls, msg):
        if cls.mac_log == True:
            pass
    
    @classmethod
    def log_net(cls, msg):
        if cls.net_log == True:
            pass
    
    @classmethod
    def log_app(cls, msg):
        if cls.app_log == True:
            pass
    
    class Rx_Log:
        def __init__(self,interface,prob,rx_sinr):
            self.interface = interface
            self.prob = prob
            self.rx_sinr = rx_sinr