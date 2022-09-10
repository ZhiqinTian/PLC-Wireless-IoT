from Network import *
from Phy import *
import matplotlib.pyplot as plt

nodes = Node.creat_nodes(8)
plc_topo = Topology()
plc_topo.set_nodes(nodes)
for node in nodes:
    plc_topo.set_load(node, 1e8)
lenth = 50
l1 = 0.3*lenth
l2 = 0.4*lenth
l3 = 0.8*lenth
plc_topo.add_line(nodes[0], nodes[2], l1)
plc_topo.add_line(nodes[3], nodes[2], l2-l1)
plc_topo.add_line(nodes[4], nodes[3], l3-l2)
plc_topo.add_line(nodes[4], nodes[1], lenth-l3)
plc_topo.add_line(nodes[5], nodes[2], 10)
plc_topo.add_line(nodes[6], nodes[3], 10)
plc_topo.add_line(nodes[7], nodes[4], 10)
#plc_topo.add_line(nodes[3], nodes[4], 150)
#plc_topo.add_line(nodes[5], nodes[1], 170)
plc_topo.show()

sm = PLC_Phy.Spectrum_Model(2e6, 30e6, 1000)
plc_channel = PLC_Channel(sm)
plc_channel.set_topo(plc_topo)
phy = PLC_Phy()
f = np.array(sm.get_subcarriers())/1e6
tx_psd = phy.psd_alloc(20)
noise_psd = np.power(
    10, (-145+53.23*np.power(sm.subcarriers*1e-6, -0.337))/10)
plt.figure(1)
plt.plot(f, 10*np.log10(noise_psd), label='noise')
plt.plot(f, 10*np.log10(tx_psd), label='tx_psd')
plt.xlabel('f MHz')
plt.ylabel('N/tx_p dBm/hz')
plt.legend()
for i in range(20):
    l1 = 0.3*(lenth+50*i)
    l2 = 0.4*(lenth+50*i)
    l3 = 0.8*(lenth+50*i)
    plc_topo.set_line_lenth(nodes[0], nodes[2], l1)
    plc_topo.set_line_lenth(nodes[3], nodes[2], l2-l1)
    plc_topo.set_line_lenth(nodes[4], nodes[3], l3-l2)
    plc_topo.set_line_lenth(nodes[4], nodes[1], lenth+50*i-l3)
    plc_channel.nodes_transfers_calc(nodes)
    transfer = plc_channel.get_transfer(nodes[0], nodes[1])
    rx_psd = tx_psd*np.square(np.power(10, transfer/10))
    snr = phy.sinr_model.snr_calc(rx_psd, noise_psd)
    plt.figure(2)
    plt.plot(f, transfer, label=str(lenth+50*i)+'m')
    plt.figure(3)
    plt.plot(f, 10*np.log10(rx_psd), label=str(lenth+50*i)+'m')
    plt.figure(4)
    plt.plot(f, 10*np.log10(snr), label=str(lenth+50*i)+'m')

nodes = Node.creat_nodes(2)
plc_topo = Topology()
plc_topo.set_nodes(nodes)
for node in nodes:
    plc_topo.set_load(node, 1e8)
lenth = 50
plc_topo.add_line(nodes[0], nodes[1], lenth)
#plc_topo.show()

sm = PLC_Phy.Spectrum_Model(2e6, 30e6, 1000)
plc_channel = PLC_Channel(sm)
plc_channel.set_topo(plc_topo)
phy = PLC_Phy()
f = np.array(sm.get_subcarriers())/1e6
tx_psd = phy.psd_alloc(20)
noise_psd = np.power(
    10, (-145+53.23*np.power(sm.subcarriers*1e-6, -0.337))/10)
plt.figure(1)
plt.plot(f, 10*np.log10(noise_psd), label='noise')
plt.plot(f, 10*np.log10(tx_psd), label='tx_psd')
plt.xlabel('f MHz')
plt.ylabel('N/tx_p dBm/hz')
plt.legend()
for i in range(20):
    plc_topo.set_line_lenth(nodes[0], nodes[1], lenth+50*i)
    plc_channel.nodes_transfers_calc(nodes)
    transfer = plc_channel.get_transfer(nodes[0], nodes[1])
    rx_psd = tx_psd*np.square(np.power(10, transfer/10))
    snr = phy.sinr_model.snr_calc(rx_psd, noise_psd)
    plt.figure(2)
    plt.plot(f, transfer, linewidth=0.5,color='black',linestyle='--')
    plt.figure(3)
    plt.plot(f, 10*np.log10(rx_psd), label=str(lenth+50*i)+'m')
    plt.figure(4)
    plt.plot(f, 10*np.log10(snr), label=str(lenth+50*i)+'m')

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
plt.show()
