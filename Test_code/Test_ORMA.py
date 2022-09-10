from ORM.ORMA_new import *
from Phy.PLC_Phy import PLC_Phy
from Phy.RF_Channel import RF_Channel
from Phy.RF_Phy import RF_Phy
from scipy.special import iv
import matplotlib.pyplot as plt


x=np.linspace(0.001,1,1000)
e_2=0.1
s=1
hf_2=1
a=ncx2.pdf(x, 2, hf_2*2/e_2, scale=e_2/2)/ncx2.cdf(x, 2, hf_2*2/e_2, scale=e_2/2)
#b=np.log(lognorm.cdf(np.exp(x), s, scale=1 / np.exp(s ** 2 / 2)))
#b=1/e_2*np.exp(-(hf_2+x)/e_2)*iv(0,2*np.sqrt(hf_2*x/e_2**2))
plt.plot(x,a)
plt.show()
exit(0)
sm=PLC_Phy.Spectrum_Model(num_bands=20)
N = sm.num_subcarriers
np.random.seed(1)
orma=ORMA(1e-3,0.1,sm)


rf_sm = RF_Phy.Spectrum_Model(2.4e9, 2e6, 22e6, 20)
rf_channel = RF_Channel(rf_sm)
hf_2=rf_channel.multipath_model.power_gain(rf_sm)
print('hf_2',np.round(hf_2,2))
s_p=np.ones(N)*1
sub_band = sm.spacing
P_w_all = 100
P_p_all = 100
# snr_w_ref = np.round(np.random.uniform(10, 35, size=N), 2)
Index=list(range(N))
snr_p_ref = np.round(np.random.uniform(40, 60, size=N), 2)
cnr_w = np.power(10, 5/ 10) * N / P_w_all
cnr_p = np.power(10, snr_p_ref / 10) * N / P_p_all

#result=orma.Outage_Rate_Max_Algorithm(cnr_w,hf_2,cnr_p,s_p,P_w_all,P_p_all)
div,mul = orma.AADS(cnr_w,hf_2,cnr_p,s_p,P_w_all,P_p_all)
#_,_,R_aads,iter_num=orma.optimal_power_allocation(cnr_w,hf_2,cnr_p,s_p,P_w_all,P_p_all, div, mul)
#_,_,R_mul,iter_num=orma.optimal_power_allocation(cnr_w,hf_2,cnr_p,s_p, P_w_all, P_p_all, [], Index)

_,_,R_div,iter_num=orma.optimal_power_allocation(cnr_w,hf_2,cnr_p,s_p,P_w_all,P_p_all, Index, [])
#print(div)
# print('opt',R_aads)
# print('mul',R_mul)
print('div',R_div)