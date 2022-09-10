from Sim.Simpy_Core import *
from Sim.Simulator import *
from Network.Packet import *
from .Nist_Error_Model import *
from .Link_Adap import *
from .MCS import *
from Sim.Time import *
import numpy as np
import math


class PLC_Noise:
    GOOD = 0
    BAD = 1


class PLC_Phy:

    def __init__(self):
        self.symbol_duration = MicroSecond(4)
        self.spectrum_model = self.Spectrum_Model()
        self.defult_mcs = MCS(BPSK, CodeRate_1_2)
        self.sinr_model = self.SINR_Model()
        self.error_rate_model = Nist_ErrorRate_Model()

    def set_mcs(self, mcs):
        self.defult_mcs = mcs

    def set_noise(self, case):
        self.sinr_model.set_noise(case, self.spectrum_model)

    def channel_est(self, tx_psd, rx_psd, interfer_psd):
        sinr_m = self.sinr_model
        sinr_m.set_interference(interfer_psd)
        rx_sinr = sinr_m.sinr_calc(rx_psd)
        h = np.sqrt(rx_psd/tx_psd)
        N = rx_psd/rx_sinr
        est_para = self.Channel_Para(h, N)
        return est_para, rx_sinr

    def get_symbols_num(self, phy_frame):
        bits_per_ofdm_symbol = np.sum(
            [MCS.get_bits_per_symbol(mcs) for mcs in phy_frame.preamble.data_mcs])
        num_symbols = phy_frame.preamble.symbols_num + \
            math.ceil(phy_frame.size.bits/bits_per_ofdm_symbol)
        return num_symbols

    def tx_duration_calc(self, phy_frame):
        num_symbols = self.get_symbols_num(phy_frame)
        duration = MilliSecond(num_symbols*self.symbol_duration.sim_time)
        return duration, num_symbols

    def robo_psd_alloc(self, tx_power):
        tx_psd = Link_Adap.pre_psd_alloc(self.spectrum_model, tx_power)
        return tx_psd

    def broadcast_mode(self, preamble, tx_power):
        tx_psd = self.robo_psd_alloc(tx_power)
        preamble.p_mcs = MCS(BPSK, CodeRate_1_2)
        preamble.data_mcs = [self.defult_mcs] * \
            self.spectrum_model.num_subcarriers
        preamble.carriers = [True]*self.spectrum_model.num_subcarriers
        return tx_psd

    def unicast_mode(self, preamble, tx_power, to_para_est):
        tx_psd = Link_Adap.pre_psd_alloc(self.spectrum_model, tx_power)
        h = to_para_est.h
        N = to_para_est.N
        sinr_est = tx_psd*np.square(h)/N
        preamble.p_mcs = MCS(BPSK, CodeRate_1_2)
        preamble.data_mcs = Link_Adap.adap_mcs_selec(sinr_est)
        preamble.carriers = [True]*self.spectrum_model.num_subcarriers
        return tx_psd

    class Spectrum_Model:
        down_bound = None
        up_bound = None
        num_subcarriers = None
        spacing = None
        subcarriers = None

        def __init__(self, down_bound=2e6, up_bound=30e6, num_bands=112):
            if not self.down_bound:
                self.set(down_bound, up_bound, num_bands)

        @ classmethod
        def set(cls, down_bound, up_bound, num_bands):
            cls.down_bound = down_bound
            cls.up_bound = up_bound
            cls.num_subcarriers = num_bands
            cls.spacing = (up_bound-down_bound)/num_bands
            cls.subcarriers = np.linspace(
                cls.down_bound+0.5*cls.spacing, cls.up_bound-0.5*cls.spacing, cls.num_subcarriers)

        @ classmethod
        def get_subcarriers(cls):
            return cls.subcarriers

    class SINR_Model:
        a = [-145, -140]
        b = [53.23, 38.75]
        c = [-0.337, -0.72]
        defult_noise_case = PLC_Noise.GOOD

        def __init__(self):
            self.noise_psd = None
            self.interference_psd = None

        def set_noise(self, case, sm):
            a = self.a[case]
            b = self.b[case]
            c = self.c[case]
            noise_psd = np.power(
                10, (a+b*np.power(sm.subcarriers*1e-6, c))/10)
            self.noise_psd = np.array([min(val, 1e-12) for val in noise_psd])
            # print('noise',noise_psd)

        def set_interference(self, interference_psd):
            self.interference_psd = interference_psd

        def snr_calc(self, rx_psd, noise_psd=None):
            if noise_psd is None:
                if self.noise_psd is None:
                    self.set_noise(self.defult_noise_case,
                                   PLC_Phy.Spectrum_Model)
                noise_psd = self.noise_psd
            snr = rx_psd/noise_psd
            return snr

        def sinr_calc(self, rx_psd, noise_psd=None):
            if noise_psd is None:
                if self.noise_psd is None:
                    self.set_noise(self.defult_noise_case,
                                   PLC_Phy.Spectrum_Model)
                noise_psd = self.noise_psd
            sinr = rx_psd/(noise_psd+self.interference_psd)
            return sinr

    class Channel_Para:
        def __init__(self, h, N):
            self.h = h
            self.N = N

    class Helper:
        def __init__(self):
            pass
