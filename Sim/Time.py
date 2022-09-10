class Second:

    def __init__(self, value):
        self.value = value
        self.type = 's'

    @property
    def sim_time(self):
        return self.value*1000

    def to_us(self):
        return MicroSecond(self.value*1000000)

    def to_ms(self):
        return MilliSecond(self.value*1000)


class MilliSecond:

    def __init__(self, value):
        self.value = value
        self.type = 'ms'

    @property
    def sim_time(self):
        return self.value

    def to_us(self):
        return self.value*1000

    def to_s(self):
        return Second(self.value/1000)


class MicroSecond:
    def __init__(self, value):
        self.value = value
        self.type = 'us'

    @property
    def sim_time(self):
        return self.value/1000

    def to_ms(self):
        return MilliSecond(self.value/1000)

    def to_s(self):
        return Second(self.value/1000000)
