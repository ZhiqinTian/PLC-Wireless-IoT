from Phy.Adap_Power_Allocat import *
import matplotlib.pyplot as plt
from Phy.PLC_Phy import *

sm = PLC_Phy.Spectrum_Model(10e6, 30e6, 100)
adap_pa = Adap_PA(m=1.2, s=1.5, p_error=1e-3, sm=sm)
N = 100
P_w_all = 100
P_p_all = 100
snr_w_ref = np.round(np.clip(np.random.normal(20, 15, size=N), 10, 30), 2)
snr_p_ref = np.round(np.clip(np.random.normal(20, 15, size=N), 10, 30), 2)

cnr_w = np.power(10, snr_w_ref / 10) * N / P_w_all
cnr_p = np.power(10, snr_p_ref / 10) * N / P_p_all
P_w1, P_p1, R_div = adap_pa.div_PA(cnr_w, P_w_all, cnr_p, P_p_all)
P_w2, P_p2, R_reuse = adap_pa.reuse_PA(cnr_w, P_w_all, cnr_p, P_p_all)
flag, p_w3,P_p3,R_adap = adap_pa.adap_PA(cnr_w, P_w_all, cnr_p, P_p_all)
#print('P_w', np.round(P_w1, 2))
#print('P_p', np.round(P_p1, 2))
print('div:', R_div*1e-6)
print('reuse:',R_reuse*1e-6)
print('adap:',R_adap*1e-6)

# plt.figure(2)
# plt.bar(np.linspace(1, N, N), snr_p_ref)
# plt.xlabel('subcarriers')
# plt.ylabel('snr_p_ref')
# plt.show()
