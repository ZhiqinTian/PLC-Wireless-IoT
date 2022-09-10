from Sim.Simpy_Core import *
from Sim.Simulator import *
from Network.Packet import *
from .Nist_Error_Model import *
from .Link_Adap import *
from .MCS import *
from Sim.Time import *
import numpy as np
import math


class RF_Noise:
    GOOD = 0
    BAD = 1


class RF_Phy:

    def __init__(self):
        self.symbol_duration = MicroSecond(4)
        self.spectrum_model = self.Spectrum_Model()
        self.defult_mcs = MCS(BPSK, CodeRate_1_2)
        self.sinr_model = self.SINR_Model()
        self.rx_sinrs = []
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
        tx_psd = Link_Adap.ave_psd_alloc(self.spectrum_model, tx_power)
        return tx_psd

    def broadcast_mode(self, preamble, tx_power):
        tx_psd = self.robo_psd_alloc(tx_power)
        preamble.p_mcs = MCS(BPSK, CodeRate_1_2)
        preamble.data_mcs = [self.defult_mcs] * \
            self.spectrum_model.num_subcarriers
        preamble.carriers = [True]*self.spectrum_model.num_subcarriers
        return tx_psd

    def unicast_mode(self, preamble, tx_power, to_para_est):
        tx_psd = Link_Adap.ave_psd_alloc(self.spectrum_model, tx_power)
        h = to_para_est.h
        N = to_para_est.N
        sinr_est = tx_psd*np.square(h)/N
        preamble.p_mcs = MCS(BPSK, CodeRate_1_2)
        preamble.data_mcs = Link_Adap.adap_mcs_selec(sinr_est)
        preamble.carriers = [True]*self.spectrum_model.num_subcarriers
        return tx_psd

    class Spectrum_Model:
        base_freq = None
        down_bound = None
        up_bound = None
        num_subcarriers = None
        spacing = None
        subcarriers = None

        def __init__(self, base_freq=2.4e9, down_bound=2e6, up_bound=30e6, num_bands=112):
            if not self.down_bound:
                self.set(base_freq, down_bound, up_bound, num_bands)

        @ classmethod
        def set(cls, base_freq, down_bound, up_bound, num_bands):
            cls.base_freq = base_freq
            cls.down_bound = base_freq+down_bound
            cls.up_bound = base_freq+up_bound
            cls.num_subcarriers = num_bands
            cls.spacing = (up_bound-down_bound)/num_bands
            cls.subcarriers = np.linspace(
                cls.down_bound+0.5*cls.spacing, cls.up_bound-0.5*cls.spacing, cls.num_subcarriers)

        @ classmethod
        def get_subcarriers(cls):
            return cls.subcarriers

    class SINR_Model:
        n = [4e-18, 1e-17]
        defult_noise_case = RF_Noise.GOOD

        def __init__(self):
            self.noise_psd = None
            self.interference_psd = None

        def set_noise(self, case, sm):
            self.noise_psd = np.array([self.n[case]]*sm.num_subcarriers)

        def set_interference(self, interference_psd):
            self.interference_psd = interference_psd

        def snr_calc(self, rx_psd, noise_psd=None):
            if noise_psd is None:
                if self.noise_psd is None:
                    self.set_noise(self.defult_noise_case,
                                   RF_Phy.Spectrum_Model)
                noise_psd = self.noise_psd
            snr = rx_psd/noise_psd
            return snr

        def sinr_calc(self, rx_psd, noise_psd=None):
            if noise_psd is None:
                if self.noise_psd is None:
                    self.set_noise(self.defult_noise_case,
                                   RF_Phy.Spectrum_Model)
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
