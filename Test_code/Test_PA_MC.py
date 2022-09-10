from ORM.ORMA import *
import matplotlib.pyplot as plt
from Network import *
from Phy import *

nodes = Node.creat_nodes(26)
topo = Topology()
pos = [
    [0, 0],
    [5, 0],
    [10, 0],
    [15, 0],
    [20, 0],
    [25, 0],
    [30, 0],
    [35, 0],
    [40, 0],
    [45, 0],

    [5, 4],
    [10, 3],
    [15, 1],
    [20, 7],
    [25, 10],
    [30, 2],
    [35, 6],
    [40, 5.8],

    [5, -4.5],
    [10, -6],
    [15, -1.7],
    [20, -5],
    [25, -2],
    [30, -3.5],
    [35, -6.9],
    [40, -3.8],
]
topo.add_nodes(nodes, pos)

topo.add_line(nodes[0], nodes[1])
topo.add_line(nodes[1], nodes[2])
topo.add_line(nodes[2], nodes[3])
topo.add_line(nodes[3], nodes[4])
topo.add_line(nodes[4], nodes[5])
topo.add_line(nodes[5], nodes[6])
topo.add_line(nodes[6], nodes[7])
topo.add_line(nodes[7], nodes[8])
topo.add_line(nodes[8], nodes[9])

topo.add_line(nodes[1], nodes[10])
topo.add_line(nodes[2], nodes[11])
topo.add_line(nodes[3], nodes[12])
topo.add_line(nodes[4], nodes[13])
topo.add_line(nodes[5], nodes[14])
topo.add_line(nodes[6], nodes[15])
topo.add_line(nodes[7], nodes[16])
topo.add_line(nodes[8], nodes[17])

topo.add_line(nodes[1], nodes[18])
topo.add_line(nodes[2], nodes[19])
topo.add_line(nodes[3], nodes[20])
topo.add_line(nodes[4], nodes[21])
topo.add_line(nodes[5], nodes[22])
topo.add_line(nodes[6], nodes[23])
topo.add_line(nodes[7], nodes[24])
topo.add_line(nodes[8], nodes[25])
topo.show()

plc_sm = PLC_Phy.Spectrum_Model(5e6, 25e6, 112)
plc_channel = PLC_Channel(plc_sm)
plc_channel.set_topo(topo)
plc_channel.nodes_transfers_calc(nodes)

rf_sm = RF_Phy.Spectrum_Model(2.4e9, 5e6, 25e6, 112)
rf_channel = RF_Channel(rf_sm)
rf_channel.set_topo(topo)

plc_phy = PLC_Phy()
rf_phy = RF_Phy()

f = np.array(plc_sm.get_subcarriers())
transfer1 = plc_channel.get_transfer(nodes[0], nodes[9])
transfer2 = rf_channel.pathloss_calc(nodes[0], nodes[9])

plc_tx_psd = Link_Adap.ave_psd_alloc(plc_sm, 112)
plc_noise_psd = np.power(
    10, (-140 + 53.23 * np.power(plc_sm.subcarriers * 1e-6, -0.337)) / 10
)
plc_rx_psd = plc_tx_psd * np.square(np.power(10, transfer1 / 10))
plc_snr = plc_phy.sinr_model.snr_calc(plc_rx_psd, plc_noise_psd)

rf_tx_psd = Link_Adap.ave_psd_alloc(rf_sm, 100)
rf_noise_psd = np.array([4e-18] * rf_sm.num_subcarriers)
rf_rx_psd = rf_tx_psd * np.square(np.power(10, transfer2 / 10))
rf_snr = rf_phy.sinr_model.snr_calc(rf_rx_psd, rf_noise_psd)
plt.figure(1)
plt.plot(f*1e-6, transfer1, label="PLC channel fading")
plt.xlabel('f /MHz',fontsize=15)
plt.ylabel('|hp| /dB',fontsize=15)
plt.grid(linestyle=':')
#plt.plot(f, transfer2, label="Wireless")
plt.legend()
plt.figure(2)
plt.plot(f*1e-6, 10 * np.log10(plc_snr), label="PLC")
# plt.xlabel('f /MHz',fontsize=15)
# plt.ylabel('SNR /dB',fontsize=15)
# plt.grid(linestyle=':')
# plt.legend()
# plt.figure(3)
plt.plot(f*1e-6, 10 * np.log10(rf_snr), color='red',label="Wireless")
plt.xlabel('f /MHz',fontsize=15)
plt.ylabel('SNR /dB',fontsize=15)
plt.grid(linestyle=':')
plt.legend()
plt.show()
