from ORM.ORMA import *
import matplotlib.pyplot as plt
from Network import *
from Phy import *
import time

np.random.seed(1)
N_room=10
N_room_node=6
topo=indoor_senerio(N_room = N_room,N_room_node = N_room_node,Level_room = 3,N_room_hor = 5,
N_room_ver = 2,L_room = 10,H_room = 15,H_room_spacing=5)
topo.show()
nodes=topo.nodes
plc_sm = PLC_Phy.Spectrum_Model(2e6, 22e6, 100)
plc_channel = PLC_Channel(plc_sm)
plc_channel.set_topo(topo)

nodes_calc=[]
for i in range(N_room):
    i1=i*N_room_node+1
    i2=(i+1)*N_room_node
    nodes_calc=nodes_calc+nodes[i1:i2]
t1=time.time()
plc_channel.nodes_transfers_calc(nodes_calc)
t2=time.time()
print('t60',t2-t1)
plc_transfers=plc_channel.transfers
np.save('plc_transfers50.npy', plc_transfers)
#plc_channel.transfers=np.load('plc_transfers.npy',allow_pickle=True).item()
rf_sm = RF_Phy.Spectrum_Model(2.4e9, 2e6, 22e6, 100)
rf_channel = RF_Channel(rf_sm)
rf_channel.set_topo(topo)

plc_phy = PLC_Phy()
rf_phy = RF_Phy()

f1 = np.array(plc_sm.get_subcarriers())
f2 = np.array(rf_sm.get_subcarriers())-2.4e9
transfer1 = plc_channel.transfer_calc(nodes[24],nodes[50])
transfer2 = rf_channel.pathloss_calc(nodes[24], nodes[50])*np.ones(rf_sm.num_subcarriers)

plc_tx_psd = Link_Adap.ave_psd_alloc(plc_sm, 100)
plc_noise_psd = np.power(
    10, (-145 + 53.23 * np.power(plc_sm.subcarriers * 1e-6, -0.337)) / 10
)
plc_rx_psd = plc_tx_psd * np.square(np.power(10, transfer1 / 10))
plc_snr = plc_phy.sinr_model.snr_calc(plc_rx_psd, plc_noise_psd)

rf_tx_psd = Link_Adap.ave_psd_alloc(rf_sm, 100)
rf_noise_psd = np.array([2e-17] * rf_sm.num_subcarriers)
rf_rx_psd = rf_tx_psd * np.square(np.power(10, transfer2 / 10))
rf_snr = 10 * np.log10(rf_phy.sinr_model.snr_calc(rf_rx_psd, rf_noise_psd))
plt.figure(1)
plt.plot(f1, transfer1, label="PLC")
plt.plot(f2, transfer2, label="Wireless")
plt.legend()
plt.figure(2)
plt.plot(f1, 10 * np.log10(plc_snr), label="PLC")
plt.plot(f2, rf_snr, label="Wireless")
#plt.plot(f2, 10*np.log10(rf_snr), label="Wireless_before")
#plt.plot(f, 10 * np.log10(rf_snr), label="RF")
# plt.legend()
plt.show()
