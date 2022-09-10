from Network import *
from Phy import *
import matplotlib.pyplot as plt

nodes = Node.creat_nodes(2)
topo = Topology()
pos = [[0, 0], [5, 0]]
topo.set_nodes(nodes, pos)
lenth = 5
topo.add_line(nodes[0], nodes[1])
topo.set_barrier(nodes[0], nodes[1], topo.Barrier_Env.GOOD)
topo.show()

plc_sm = PLC_Phy.Spectrum_Model(2e6, 30e6, 110)
plc_channel = PLC_Channel(plc_sm)
plc_channel.set_topo(topo)
plc_phy = PLC_Phy()

rf_sm = RF_Phy.Spectrum_Model(2.4e9, 2e6, 30e6, 110)
rf_channel = RF_Channel(rf_sm)
rf_channel.set_topo(topo)
rf_phy = RF_Phy()

f = np.array(plc_sm.get_subcarriers())/1e6
plc_tx_psd = Link_Adap.ave_psd_alloc(plc_sm, 100)
plc_noise_psd = np.power(
    10, (-145+53.23*np.power(plc_sm.subcarriers*1e-6, -0.337))/10)
# plc_noise_psd = np.power(
#     10, (-140+38.75*np.power(plc_sm.subcarriers*1e-6, -0.72))/10)
#plc_noise_psd = np.array([min(val, 1e-12) for val in plc_noise_psd])
plt.figure(1)
plt.plot(f, 10*np.log10(plc_noise_psd), label='plc_noise')
plt.plot(f, 10*np.log10(plc_tx_psd), label='plc_tx_psd')
plt.xlabel('f MHz')
plt.ylabel('N/tx_p dBm/hz')
plt.legend()

rf_tx_psd = Link_Adap.ave_psd_alloc(rf_sm, 100)
rf_noise_psd = np.array([4e-18]*rf_sm.num_subcarriers)
plt.figure(5)
plt.plot(f, 10*np.log10(rf_noise_psd), label='plc_noise')
plt.plot(f, 10*np.log10(rf_tx_psd), label='plc_tx_psd')
plt.xlabel('f MHz')
plt.ylabel('N/tx_p dBm/hz')
plt.legend()

for i in range(20):
    topo.set_node_pos(nodes[1], lenth+2*i)
    topo.set_line_lenth(nodes[0], nodes[1], lenth+2*i)
    plc_channel.nodes_transfers_calc(nodes)
    transfer1 = plc_channel.get_transfer(nodes[0], nodes[1])
    # if i == 19:
    #     print(f[0], transfer1[0], f[-1], transfer1[-1])
    #     a1 =(transfer1[-1]-transfer1[0])/1e6/(-868.59)/(f[-1]-f[0])
    #     a0 =(transfer1[0]/(-868.59)-a1*f[0])
    #     print(a1,a0)
    plc_rx_psd = plc_tx_psd*np.square(np.power(10, transfer1/10))
    plc_snr = plc_phy.sinr_model.snr_calc(plc_rx_psd, plc_noise_psd)

    transfer2 = rf_channel.pathloss_calc(nodes[0], nodes[1])
    rf_rx_psd = rf_tx_psd*np.square(np.power(10, transfer2/10))
    rf_snr = rf_phy.sinr_model.snr_calc(rf_rx_psd, rf_noise_psd)

    plt.figure(2)
    plt.plot(f, transfer1, label=str(lenth+2*i)+'m')
    plt.figure(3)
    plt.plot(f, 10*np.log10(plc_rx_psd), label=str(lenth+2*i)+'m')
    plt.figure(4)
    plt.plot(f, 10*np.log10(plc_snr), label=str(lenth+2*i)+'m')

    plt.figure(6)
    plt.plot(f, transfer2, label=str(lenth+2*i)+'m')
    plt.figure(7)
    plt.plot(f, 10*np.log10(rf_rx_psd), label=str(lenth+2*i)+'m')
    plt.figure(8)
    plt.plot(f, 10*np.log10(rf_snr), label=str(lenth+2*i)+'m')

plt.figure(2)
plt.xlabel('f MHz')
plt.ylabel('H dB')
plt.legend()
plt.figure(3)
plt.xlabel('f MHz')
plt.ylabel('rx_psd dBm/Hz')
plt.legend()
plt.figure(4)
plt.plot(f, 4*np.ones(len(f)), linewidth=2, linestyle='-.')
plt.xlabel('f MHz')
plt.ylabel('snr dB')
plt.legend()

plt.figure(6)
plt.xlabel('f MHz')
plt.ylabel('H dB')
plt.legend()
plt.figure(7)
plt.xlabel('f MHz')
plt.ylabel('rx_psd dBm/Hz')
plt.legend()
plt.figure(8)
plt.plot(f, 4*np.ones(len(f)), linewidth=2, linestyle='-.')
plt.xlabel('f MHz')
plt.ylabel('snr dB')
plt.legend()
plt.show()
