import os
import numpy as np
import matlab.engine
from Link.PLC_Interface import PLC_Interfaces
import time

class PLC_Channel:
    def __init__(self, sm=None):
        self.plc_topo = None
        self.spectrum_model = sm
        self.transfers = {}
        self.plc_mat_api = PLC_MATLAB_API()

    def set_topo(self, topo):
        self.plc_topo = topo

    def get_transfer(self, tx_node, rx_node):
        return self.transfers[tx_node.id][rx_node.id]

    def transmission_calc(self, tx_interface, rx_interface, tx_psd):
        transfer_vector = self.get_transfer(
            tx_interface.node, rx_interface.node)
        rx_psd = tx_psd*np.square(np.power(10, transfer_vector/10))
        return rx_psd

    def transfer_calc(self, tx_node, rx_node):
        self.plc_topo.bone_paths_refresh()
        n, units_struct = self.plc_topo.units_struct_generate(tx_node, rx_node)
        f = self.spectrum_model.get_subcarriers()
        return self.plc_mat_api.transfer_calc(f, n, units_struct)

    def nodes_transfers_calc(self, nodes):
        self.transfers = {}
        self.plc_topo.bone_paths_refresh()
        num=0
        for tx in nodes:
            for rx in nodes:
                if tx != rx:
                    num+=1
                    t1=time.time()
                    n, units_struct = self.plc_topo.units_struct_generate(
                        tx, rx)
                    f = self.spectrum_model.get_subcarriers()
                    transfer = self.plc_mat_api.transfer_calc(
                        f, n, units_struct)
                    if tx.id not in self.transfers:
                        self.transfers[tx.id] = {}
                    self.transfers[tx.id][rx.id] = transfer
                    t2=time.time()
                    print(num,':',t2-t1)

    def interfaces_transfers_calc(self):
        self.transfers = {}
        plc_interfaces = PLC_Interfaces.get_all()
        self.plc_topo.bone_paths_refresh()
        for tx in plc_interfaces:
            for rx in plc_interfaces:
                if tx != rx:
                    n, units_struct = self.plc_topo.units_struct_generate(
                        tx.node, rx.node)
                    f = self.spectrum_model.get_subcarriers()
                    transfer = self.plc_mat_api.transfer_calc(
                        f, n, units_struct)
                    if tx.node.id not in self.transfers:
                        self.transfers[tx.node.id] = {}
                    self.transfers[tx.node.id][rx.node.id] = transfer

    def transmission_range_calc(self, tx_interface, tx_psd, threshold):
        trans_range = []
        interfaces = PLC_Interfaces.get_all()
        for rx_interface in interfaces:
            if rx_interface != tx_interface:
                rx_psd = self.transmission_calc(
                    tx_interface, rx_interface, tx_psd)
                spacing = self.spectrum_model.spacing
                rx_power = spacing*np.sum(rx_psd)
                if rx_power >= threshold:
                    trans_range.append(rx_interface)
        return trans_range

    def channel_calc(self):
        self.interfaces_transfers_calc()
        interfaces = PLC_Interfaces.get_all()
        for interface in interfaces:
            interface.update_cca_neighbors()


class PLC_MATLAB_API:
    def __init__(self):
        self.current_path = os.getcwd()
        self.eng = matlab.engine.start_matlab()
        self.eng.cd(self.current_path+'/Phy/plc_channel_matlab')

    def data_format(self, units_struct):

        def to_mat(unit, level_num):
            def zero_padding(row, lenth):
                row += [0]*(lenth-len(row))

            max_level_line_num = max(unit['level_cable_num'])
            for level in range(level_num):
                zero_padding(unit['level_load_type']
                             [level], max_level_line_num)
                zero_padding(unit['level_cable_type']
                             [level], max_level_line_num)
                zero_padding(unit['level_len'][level], max_level_line_num)
                zero_padding(unit['next_connect'][level], max_level_line_num)
                zero_padding(unit['level_node_load']
                             [level], max_level_line_num)
            unit['level_cable_num'] = matlab.double(unit['level_cable_num'])
            unit['level_load_type'] = matlab.double(unit['level_load_type'])
            unit['level_cable_type'] = matlab.double(unit['level_cable_type'])
            unit['level_len'] = matlab.double(unit['level_len'])
            unit['next_connect'] = matlab.double(unit['next_connect'])
            unit['level_node_load'] = matlab.double(unit['level_node_load'])

        for unit in units_struct:
            unit['bone1_len'] = matlab.double([unit['bone1_len']])
            unit['bone2_len'] = matlab.double([unit['bone2_len']])
            level_num = len(unit['level_cable_num'])
            if level_num:
                to_mat(unit, level_num)

    def transfer_calc(self, subcarriers, n, units_structure):
        f = matlab.double(subcarriers.tolist())
        self.data_format(units_structure)
        transfer = self.eng.transfer_calc(f, n, units_structure)
        if isinstance(transfer, float):
            return np.array([transfer])
        return np.array(transfer[0])
