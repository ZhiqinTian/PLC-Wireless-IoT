import numpy as np
from scipy.stats import nakagami,norm,gamma
from Link.RF_Interface import RF_Interfaces
from scipy.signal import savgol_filter

class RF_Channel:

    def __init__(self, sm=None):
        self.topo = None
        self.spectrum_model = sm
        self.expect_dis_range = 250
        self.pathloss_model = self.Pathloss_Model()
        self.multipath_model = self.Nakagami_m_Model()

    def set_topo(self, topo):
        self.topo = topo

    def pathloss_calc(self, tx_node, rx_node):
        propagation = self.topo.get_propagation(tx_node, rx_node)
        dis = self.topo.get_nodes_dis(tx_node, rx_node)
        pathloss = self.pathloss_model.pathloss_calc(
            dis, self.spectrum_model, propagation)
        #pathloss = savgol_filter(pathloss, 31, 5, mode= 'nearest')
        #pathloss = savgol_filter(pathloss, 31, 5, mode= 'nearest')
        return -pathloss/2

    def multipath_loss_calc(self, tx_node, rx_node):
        dis = self.topo.get_nodes_dis(tx_node, rx_node)
        multipath_transfer = self.multipath_model.transfer_calc(
            self.spectrum_model)
        return 10*np.log10(multipath_transfer)

    def transmission_calc(self, tx_interface, rx_interface, tx_psd):
        path_transfer = self.pathloss_calc(
            tx_interface.node, rx_interface.node)
        multipath_transfer = self.multipath_loss_calc(
            tx_interface.node, rx_interface.node)
        transfer = path_transfer + multipath_transfer
        rx_psd = tx_psd*np.square(np.power(10, transfer/10))
        return rx_psd, transfer

    def est_trans_calc(self, tx_interface, rx_interface, tx_psd):
        path_transfer = self.pathloss_calc(
            tx_interface.node, rx_interface.node)
        rx_psd = tx_psd*np.square(np.power(10, path_transfer/10))
        return rx_psd, path_transfer

    def transmission_range_calc(self, tx_interface, tx_psd, threshold):
        trans_range = []
        interfaces = RF_Interfaces.get_all()
        for rx_interface in interfaces:
            if rx_interface != tx_interface:
                if self.topo.get_nodes_dis(rx_interface.node, tx_interface.node) < self.expect_dis_range:
                    rx_psd, _ = self.est_trans_calc(
                        tx_interface, rx_interface, tx_psd)
                    spacing = self.spectrum_model.spacing
                    rx_power = spacing*np.sum(rx_psd)
                    if rx_power >= threshold:
                        trans_range.append(rx_interface)
        return trans_range

    def channel_calc(self):
        interfaces = RF_Interfaces.get_all()
        for interface in interfaces:
            interface.update_cca_neighbors()

    class Pathloss_Model:
        def __init__(self):
            self.n_nLOS = 3
            self.n_NLOS = 3.5
            self.Gt = 1
            self.Gr = 1
            self.d0 = 1
            self.shadow_loss = 5
            self.nLOS_loss = 5
            self.NLOS_loss = 10
            

        def pathloss_calc(self, dis, sm, propagation):
            f = sm.base_freq
            if propagation==0:
                n=self.n_nLOS
            else:
                n=self.n_NLOS
            PL0 = 10*np.log10(np.square(4*3.14*f*self.d0/3e8) /
                              (self.Gt*self.Gr))+propagation*self.NLOS_loss+(1-propagation)*self.nLOS_loss
            S = np.random.normal( 0, self.shadow_loss)
            pathloss = PL0 + 10*n * np.log10(dis/self.d0) + S
            return pathloss
            

    class Nakagami_m_Model:
        def __init__(self):
            self.oumiga = 1
            self.m = 1

        def transfer_calc(self, sm, m=None):
            if m is None:
                m = self.m
            transfer = np.array(nakagami.rvs(
                m, loc=0, scale=self.oumiga, size=sm.num_subcarriers))
            return transfer
        
        def power_gain(self, sm, m=None):
            if m is None:
                m = self.m
            gain = np.array(gamma.rvs( m, scale=1 / m, size=sm.num_subcarriers))
            return gain