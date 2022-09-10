from networkx.classes.function import neighbors
from Link.PLC_Interface import PLC_Interface, PLC_Interfaces
from Network import *
from Phy import *
import matplotlib.pyplot as plt

nodes = Node.creat_nodes(9)
plc_topo = Topology()
plc_topo.set_nodes(nodes)
for node in nodes:
    plc_topo.set_load(node, 1e8)

plc_topo.add_line(nodes[0], nodes[7], 80)
plc_topo.add_line(nodes[1], nodes[7], 100)
plc_topo.add_line(nodes[1], nodes[2], 100)
plc_topo.add_line(nodes[3], nodes[1], 150)
plc_topo.add_line(nodes[5], nodes[1], 170)
plc_topo.add_line(nodes[6], nodes[5], 200)
plc_topo.add_line(nodes[3], nodes[4], 80)
plc_topo.add_line(nodes[5], nodes[8], 100)
plc_topo.show()

sm = PLC_Phy.Spectrum_Model(0, 100e6, 2000)
plc_channel = PLC_Channel(sm)
plc_channel.set_topo(plc_topo)
for i in range(len(nodes)):
    phy = PLC_Phy()
    phy.set_channel(plc_channel)
    # phy.set_mcs(MCS.QPSK_1_2)
    phy.set_mcs(MCS.QAM16_1_2)

    interface = PLC_Interface()
    interface.set_node(nodes[i])
    interface.set_phy(phy)
    interface.set_tx_power(20)

plc_channel.channel_calc()

mac_f = MAC_Frame()
tx = PLC_Interfaces.get(nodes[0])
rx = PLC_Interfaces.get(nodes[8])
Simulator.schedule_now(tx.send, mac_f, 0)
Simulator.run()
rx_psd = plc_channel.transmission_calc(tx, rx, tx.tx_psd)
spacing = sm.spacing
rx_power = spacing*sum(rx_psd)
print('rx_power', rx_power)
neighbors = [i.node.id for i in tx.cca_neighbours]
print(neighbors)

plt.show()
