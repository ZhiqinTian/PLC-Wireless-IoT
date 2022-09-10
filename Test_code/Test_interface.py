from Link.PLC_Interface import PLC_Interface
from Network import *
from Phy import *
from Link.Arq_MAC import Arq_MAC
import matplotlib.pyplot as plt

nodes = Node.creat_nodes(9)
plc_topo = Topology()
plc_topo.set_nodes(nodes)
for node in nodes:
    plc_topo.set_load(node, 1e8)

plc_topo.add_line(nodes[0], nodes[7], 200)
plc_topo.add_line(nodes[1], nodes[7], 80)
plc_topo.add_line(nodes[1], nodes[2], 100)
plc_topo.add_line(nodes[3], nodes[1], 150)
plc_topo.add_line(nodes[5], nodes[1], 170)
plc_topo.add_line(nodes[6], nodes[5], 200)
plc_topo.add_line(nodes[3], nodes[4], 80)
plc_topo.add_line(nodes[5], nodes[8], 200)
plc_topo.show()

sm = PLC_Phy.Spectrum_Model(2e6, 30e6, 300)
plc_channel = PLC_Channel(sm)
plc_channel.set_topo(plc_topo)

tx_phy = PLC_Phy()
tx_phy.set_channel(plc_channel)
tx_phy.set_mcs(MCS(BPSK, CodeRate_1_2))
# tx_phy.set_mcs(MCS.QAM16_1_2)

rx_phy = PLC_Phy()
rx_phy.set_channel(plc_channel)
rx_phy.set_mcs(MCS(BPSK, CodeRate_1_2))
# rx_phy.set_mcs(MCS.QAM16_1_2)

tx = PLC_Interface()
tx.set_node(nodes[0])
tx.set_phy(tx_phy)
tx.set_tx_power(20)

rx = PLC_Interface()
rx.set_node(nodes[7])
rx.set_phy(rx_phy)
rx.set_tx_power(20)

plc_channel.channel_calc()

mac_f = MAC_Frame()
tx_mode = PLC_Interface.Tx_Mode.Broadcast
tx_para = Arq_MAC.Tx_Para(tx_mode)
Simulator.schedule_now(tx.send, mac_f, tx_para)
Simulator.run()

transfer = plc_channel.get_transfer(nodes[0], nodes[7])
# print(tx.cca_neighbours)
f = np.array(sm.get_subcarriers())
plt.figure()
plt.plot(f, transfer)
plt.show()
