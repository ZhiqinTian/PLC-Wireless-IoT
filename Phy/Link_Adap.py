from .MCS import *
import numpy as np


class Link_Adap:
    mcs_table = {0: (BPSK, CodeRate_1_2), 1: (BPSK, CodeRate_2_3), 2: (BPSK, CodeRate_3_4),
                 3: (QPSK, CodeRate_1_2), 4: (QPSK, CodeRate_2_3), 5: (QPSK, CodeRate_3_4),
                 6: (QAM8, CodeRate_2_3), 7: (QAM8, CodeRate_3_4),
                 8: (QAM16, CodeRate_2_3), 9: (QAM16, CodeRate_3_4),
                 10: (QAM32, CodeRate_2_3), 11: (QAM32, CodeRate_3_4),
                 12: (QAM64, CodeRate_2_3), 13: (QAM64, CodeRate_3_4),
                 14: (QAM128, CodeRate_2_3), 15: (QAM128, CodeRate_3_4),
                 16: (QAM256, CodeRate_2_3), 17: (QAM256, CodeRate_3_4), }
    sinr_threshold = np.array([3.9, 5.7, 6.8,
                               6.9, 8.7, 9.8,
                               12.2, 13.4,
                               15.4, 16.6,
                               18.5, 19.7,
                               21.3, 22.6,
                               24.2, 25.5,
                               27.1, 28.5])

    @ classmethod
    def ave_psd_alloc(cls, sm, tx_power):
        # by average
        num_subcarriers = sm.num_subcarriers
        spacing = sm.spacing
        tx_psd = np.array(
            [tx_power/spacing/num_subcarriers] * num_subcarriers)
        return tx_psd

    @classmethod
    def pre_psd_alloc(cls, sm, tx_power):
        a0 = 6.259e-4
        a1 = 1.397e-9
        def h_pre(f, lenth): return np.exp(-(a0+a1*f)*lenth)

        spacing = sm.spacing
        exp_dis = 100
        factor = 1/np.square(h_pre(sm.get_subcarriers(), exp_dis))
        tx_psd = tx_power*(factor/np.sum(factor))/spacing
        return tx_psd

    @ classmethod
    def adap_psd_alloc(cls):
        pass

    @ classmethod
    def adap_mcs_selec(cls, sinr=None):
        mcs = []
        if sinr:
            for sub_sinr in 10*np.log10(sinr):
                select_set = np.where(sub_sinr > cls.sinr_threshold)
                best_index = 0
                if select_set:
                    best_index = np.max(select_set)
                sub_mcs = MCS(*cls.mcs_table[best_index])
                mcs.append(sub_mcs)
        return mcs
