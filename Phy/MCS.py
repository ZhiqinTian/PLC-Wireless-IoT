BPSK = 'BPSK'
QPSK = 'QPSK'
QAM8 = 'QAM8'
QAM16 = 'QAM16'
QAM32 = 'QAM32'
QAM64 = 'QAM64'
QAM128 = 'QAM128'
QAM256 = 'QAM256'

CodeRate_1_2 = '1_2'
CodeRate_2_3 = '2_3'
CodeRate_3_4 = '3_4'


class MCS:
    def __init__(self, modulation, code_rate):
        self.modulation = modulation
        self.code_rate = code_rate

    @classmethod
    def get_code_rate(cls, mcs):
        code_rate = None
        if mcs.code_rate == CodeRate_1_2:
            code_rate = 1/2
        elif mcs.code_rate == CodeRate_2_3:
            code_rate = 2/3
        elif mcs.code_rate == CodeRate_3_4:
            code_rate = 3/4
        return code_rate

    @classmethod
    def get_bits_per_symbol(cls, mcs):
        bits = None
        if mcs.modulation == BPSK:
            bits = 1
        elif mcs.modulation == QPSK:
            bits = 2
        elif mcs.modulation == QAM8:
            bits = 3
        elif mcs.modulation == QAM16:
            bits = 4
        elif mcs.modulation == QAM32:
            bits = 5
        elif mcs.modulation == QAM64:
            bits = 6
        elif mcs.modulation == QAM128:
            bits = 7
        elif mcs.modulation == QAM256:
            bits = 8
        return bits

    @classmethod
    def get_Ms(cls):
        return [BPSK, QPSK, QAM8, QAM16, QAM32, QAM64, QAM128, QAM256]

    @classmethod
    def get_Cs(cls):
        return [CodeRate_1_2, CodeRate_2_3, CodeRate_3_4]
