from .Time import *
from Network.Packet import *
from .Simulator import *


class Point_to_Point:

    # Constructor
    def __init__(self, sender=None, receiver=None):
        self.sender = sender
        self.receiver = receiver
        self.packet = None
        self.start = Second(0)
        self.times = 1
        self.interval = Second(1)
        Simulator.add_app(self)

    # Set app packet
    def set_packet(self, packet):
        self.packet = packet

    # set times and interval
    def set_time(self, times, interval):
        self.times = times
        self.interval = interval

    # send app packet
    def send(self):
        self.sender.send_to(self.receiver, self.packet)

    # app running
    def run(self):
        for i in range(self.times):
            Simulator.schedule(self.start.sim_time+i *
                               self.interval.sim_time, self.send())
