from Link.PLC_Interface import PLC_Interface
from Link.RF_Interface import RF_Interface
from Network import *
from Phy import *
import matplotlib.pyplot as plt

nodes = Node.creat_nodes(2)
topo = Topology()
lenth = 5
pos = [[0, 0], [lenth, 0]]
topo.set_nodes(nodes, pos)
topo.add_line(nodes[0], nodes[1])
topo.set_barrier(nodes[0], nodes[1], topo.Barrier_Env.GOOD)
topo.show()

plc_sm = PLC_Phy.Spectrum_Model(2e6, 30e6, 110)
plc_channel = PLC_Channel(plc_sm)
plc_channel.set_topo(topo)

plc_tx_phy = PLC_Phy()
# tx_phy.set_mcs(MCS.QPSK_1_2)
plc_tx_phy.set_mcs(MCS(BPSK, CodeRate_1_2))
plc_tx = PLC_Interface()
plc_tx.set_channel(plc_channel)
plc_tx.set_node(nodes[0])
plc_tx.set_phy(plc_tx_phy)
plc_tx.set_tx_power(20)

plc_rx_phy = PLC_Phy()
# rx_phy.set_mcs(MCS.QPSK_1_2)
plc_rx_phy.set_mcs(MCS(QPSK, CodeRate_1_2))
plc_rx = PLC_Interface()
plc_rx.set_channel(plc_channel)
plc_rx.set_node(nodes[1])
plc_rx.set_phy(plc_rx_phy)
plc_rx.set_tx_power(20)

plc_channel.channel_calc()

rf_sm = RF_Phy.Spectrum_Model(2.4e9, 2e6, 30e6, 110)
rf_channel = RF_Channel(rf_sm)
rf_channel.set_topo(topo)

rf_tx_phy = RF_Phy()
# tx_phy.set_mcs(MCS.QPSK_1_2)
rf_tx_phy.set_mcs(MCS(BPSK, CodeRate_1_2))
rf_tx = RF_Interface()
rf_tx.set_channel(rf_channel)
rf_tx.set_node(nodes[0])
rf_tx.set_phy(rf_tx_phy)
rf_tx.set_tx_power(100)

rf_rx_phy = RF_Phy()
# rx_phy.set_mcs(MCS.QPSK_1_2)
rf_rx_phy.set_mcs(MCS(BPSK, CodeRate_1_2))
rf_rx = RF_Interface()
rf_rx.set_channel(rf_channel)
rf_rx.set_node(nodes[1])
rf_rx.set_phy(rf_rx_phy)
rf_rx.set_tx_power(100)

rf_channel.channel_calc()

mac_f1 = MAC_Frame(payload_size=Byte(1000))
plc_tx_mode = plc_tx.Tx_Mode.Broadcast
plc_tx_para = plc_tx.Tx_Para(plc_tx_mode)
Simulator.schedule_now(plc_tx.send, mac_f1, plc_tx_para)

mac_f2 = MAC_Frame(payload_size=Byte(1000))
rf_tx_mode = rf_tx.Tx_Mode.Broadcast
rf_tx_para = rf_tx.Tx_Para(rf_tx_mode)
Simulator.schedule_now(rf_tx.send, mac_f2, rf_tx_para)

for i in range(100):
    time = Second(1)
    Simulator.schedule(
        time.sim_time*(i+1), topo.set_node_pos, nodes[1], [lenth+i+1, 0])
    Simulator.schedule(
        time.sim_time*(i+1), topo.set_line_lenth, nodes[0], nodes[1], lenth+i+1)
    Simulator.schedule(
        time.sim_time*(i+1)+1, plc_channel.channel_calc)
    Simulator.schedule(
        time.sim_time*(i+1)+1, rf_channel.channel_calc)
    Simulator.schedule(
        time.sim_time*(i+1)+2, plc_tx.send, mac_f1, plc_tx_para)
    Simulator.schedule(
        time.sim_time*(i+1)+2, rf_tx.send, mac_f2, rf_tx_para)

Simulator.run()
plc_prob = [i.prob for i in plc_rx.rx_log]
plc_sinr = np.array([i.rx_sinr[109] for i in plc_rx.rx_log])

rf_prob = [i.prob for i in rf_rx.rx_log]
rf_sinr = np.array([i.rx_sinr[109] for i in rf_rx.rx_log])

l = np.linspace(lenth, 100+lenth, 101)
#transfer = plc_channel.get_transfer(nodes[0], nodes[1])
#f = np.array(sm.get_subcarriers())
plt.figure(1)
plt.xlabel('line lenth /m')
plt.ylabel('sinr /dB')
plt.plot(l, 10*np.log10(plc_sinr))
plt.figure(2)
plt.ylim(0, 1.05)
plt.xlabel('line lenth /m')
plt.ylabel('rx success prob')
plt.plot(l, plc_prob)
plt.figure(3)
plt.xlabel('sinr /dB')
plt.ylabel('rx success prob')
plt.plot(10*np.log10(plc_sinr), plc_prob)

plt.figure(4)
plt.xlabel('line lenth /m')
plt.ylabel('rf_sinr /dB')
plt.plot(l, 10*np.log10(rf_sinr))
plt.figure(5)
plt.ylim(0, 1.05)
plt.xlabel('line lenth /m')
plt.ylabel('rf_rx success prob')
plt.plot(l, rf_prob)
plt.figure(6)
plt.xlabel('sinr /dB')
plt.ylabel('rf_rx success prob')
ber = np.array([rf_rx_phy.error_rate_model.ber_calc(MCS(BPSK, CodeRate_1_2), sinr) for sinr in rf_sinr])
plt.plot(np.linspace(0,1000,len(rf_sinr)), ber)
plt.show()
