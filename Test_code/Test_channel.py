from Network import *
from Phy import *
import matplotlib.pyplot as plt

nodes = Node.creat_nodes(7)
topo = Topology()
pos = [[0, 0], [6, 0], [8, 0], [15, 0], [3, 2], [8, 3], [3, -1]]
topo.set_nodes(nodes, pos)

topo.add_line(nodes[0], nodes[1])
topo.add_line(nodes[1], nodes[2])
topo.add_line(nodes[2], nodes[3])
topo.add_line(nodes[1], nodes[4])
topo.add_line(nodes[2], nodes[5])
topo.add_line(nodes[1], nodes[6])
topo.show()

plc_sm = PLC_Phy.Spectrum_Model(2e6, 30e6, 110)
plc_channel = PLC_Channel(plc_sm)
plc_channel.set_topo(topo)
plc_channel.nodes_transfers_calc(nodes)

rf_sm = RF_Phy.Spectrum_Model(2.4e9, 2e6, 30e6, 110)
rf_channel = RF_Channel(rf_sm)
rf_channel.set_topo(topo)

f = np.array(plc_sm.get_subcarriers())
transfer1 = plc_channel.get_transfer(nodes[0], nodes[1])
transfer2 = rf_channel.pathloss_calc(nodes[0], nodes[1])
plt.figure()
plt.plot(f, transfer1, label='PLC')
plt.plot(f, transfer2, label='RF')
plt.legend()
plt.show()
