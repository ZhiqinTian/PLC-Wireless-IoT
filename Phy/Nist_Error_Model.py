from .MCS import *
import numpy as np
import math


class Nist_ErrorRate_Model:
    a = {'BPSK': 1.0, 'QPSK': 1.0, 'QAM8': 0.862, 'QAM16': 0.75,
         'QAM32': 0.66, 'QAM64': 7/12, 'QAM128': 0.521, 'QAM256': 15/32}
    b = {'BPSK': 1.0, 'QPSK': 2.0, 'QAM8': 14/3, 'QAM16': 10.0,
         'QAM32': 63/3, 'QAM64': 42.0, 'QAM128': 254/3, 'QAM256': 170.0}

    @ classmethod
    def packet_success_rate_calc(cls, phy_frame, rx_sinrs):
        num_symbols = len(rx_sinrs)
        p_mcs = phy_frame.preamble.p_mcs
        data_mcs = phy_frame.preamble.data_mcs
        carriers = phy_frame.preamble.carriers
        rate = 1
        for i in range(num_symbols):
            bler = 0
            sinr = rx_sinrs[i]
            if i < phy_frame.preamble.symbols_num:
                for sub in range(len(carriers)):
                    if carriers[sub]:
                        ber = cls.ber_calc(p_mcs, sinr[sub])
                        bits = MCS.get_bits_per_symbol(p_mcs)
                        bler = 1-pow(1-ber, bits)*(1-bler)
            else:
                for sub in range(len(carriers)):
                    if carriers[sub]:
                        ber = cls.ber_calc(data_mcs[sub], sinr[sub])
                        bits = MCS.get_bits_per_symbol(data_mcs[sub])
                        bler = 1-pow(1-ber, bits)*(1-bler)
            rate *= 1-bler
        return rate

    @ classmethod
    def ber_calc(cls, mcs, sinr):
        a = cls.a[mcs.modulation]
        b = cls.b[mcs.modulation]
        ber = 0.5*a*math.erfc(math.sqrt(sinr/b))
        D = np.sqrt(4.0*ber*(1.0-ber))
        if mcs.code_rate == CodeRate_1_2:
            pe = 0.5 * (36.0 * pow(D, 10) + 211.0 * pow(D, 12) + 1404.0 * pow(D, 14) + 11633.0 * pow(D, 16) + 77433.0 * pow(
                D, 18) + 502690.0 * pow(D, 20) + 3322763.0 * pow(D, 22) + 21292910.0 * pow(D, 24) + 134365911.0 * pow(D, 26))
        elif mcs.code_rate == CodeRate_2_3:
            pe = 0.25 * (3.0 * pow(D, 6) + 70.0 * pow(D, 7) + 285.0 * pow(D, 8) + 1276.0 * pow(D, 9) + 6160.0 * pow(D, 10) +
                         27128.0 * pow(D, 11) + 117019.0 * pow(D, 12) + 498860.0 * pow(D, 13) + 2103891.0 * pow(D, 14) + 8784123.0 * pow(D, 15))
        elif mcs.code_rate == CodeRate_3_4:
            pe = 1.0 / 6.0 * (42.0 * pow(D, 5) + 201.0 * pow(D, 6) + 1492.0 * pow(D, 7) + 10469.0 * pow(D, 8) + 62935.0 * pow(
                D, 9) + 379644.0 * pow(D, 10) + 2253373.0 * pow(D, 11) + 13073811.0 * pow(D, 12) + 75152755.0 * pow(D, 13) + 428005675.0 * pow(D, 14))
        pe = min(pe, 1)
        return pe
