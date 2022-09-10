from ORM.ORMA_new import *
from Phy.PLC_Phy import PLC_Phy
from Phy.RF_Channel import RF_Channel
from Phy.RF_Phy import RF_Phy
from scipy.special import iv
import matplotlib.pyplot as plt
import time

# x=np.linspace(-2,1,1000)
# e_2=0.1
# s=1
# hf_2=1

# a=np.log(ncx2.cdf(np.exp(np.sqrt(x)), 2, hf_2*2/e_2, scale=e_2/2))
# # b=-np.log(np.exp(-3.5*x)+1)
# #b=np.log(lognorm.cdf(np.exp(x), s, scale=1 / np.exp(s ** 2 / 2)))
# #b=1/e_2*np.exp(-(hf_2+x)/e_2)*iv(0,2*np.sqrt(hf_2*x/e_2**2))
# plt.plot(x,a)
# # plt.plot(x,b)
# plt.show()
sm=PLC_Phy.Spectrum_Model(num_bands=100)
N = sm.num_subcarriers
np.random.seed(1)
orma=ORMA(1e-3,0.1,sm)

rf_sm = RF_Phy.Spectrum_Model(2.4e9, 2e6, 22e6, 100)
rf_channel = RF_Channel(rf_sm)
hf_2=rf_channel.multipath_model.power_gain(rf_sm)
hf_2=np.clip(hf_2,a_min=0.001,a_max=None)
#print('hf_2',np.round(hf_2,2))
s_p=np.ones(N)*1
sub_band = sm.spacing
P_w_all = 100
P_p_all = 100
# snr_w_ref = np.round(np.random.uniform(10, 35, size=N), 2)
Index=list(range(N))
snr_p_ref = np.round(np.random.uniform(10, 45, size=N), 2)
cnr_w = np.power(10, 50/ 10) * N / P_w_all
cnr_p = np.power(10, snr_p_ref / 10) * N / P_p_all
t1=time.time()
#result=orma.Outage_Rate_Max_Algorithm(cnr_w,hf_2,cnr_p,s_p,P_w_all,P_p_all)
div,mul = orma.AADS(cnr_w,hf_2,cnr_p,s_p,P_w_all,P_p_all)
print('div',div)
print('mul',mul)
_,_,R_aads,iter_num=orma.optimal_power_allocation(cnr_w,hf_2,cnr_p,s_p,P_w_all,P_p_all, div, mul)
_,_,R_mul,iter_num=orma.optimal_power_allocation(cnr_w,hf_2,cnr_p,s_p, P_w_all, P_p_all, [], Index)
p_w,p_p=orma.uniform_power_allocation(100,100)
p_ww,p_pw=orma.water_filling_allocation(cnr_w,hf_2,cnr_p,P_w_all,P_p_all)
p_ws,p_ps=orma.suboptimal_power_allocation(cnr_w,hf_2,cnr_p,s_p,P_w_all,P_p_all, div, mul)
R_uf=orma.outage_rate(cnr_w, hf_2, cnr_p, s_p, p_w, p_p, div, mul)
R_sub=orma.outage_rate(cnr_w, hf_2, cnr_p, s_p, p_ws, p_ps, div, mul)
_,_,R_div,iter_num=orma.optimal_power_allocation(cnr_w,hf_2,cnr_p,s_p,P_w_all,P_p_all, Index, [])
R_div_w=orma.outage_rate(cnr_w, hf_2, cnr_p, s_p, p_ww, p_pw,  Index, [])
R_mul_w=orma.outage_rate(cnr_w,hf_2,cnr_p,s_p, p_ww, p_pw, [], Index)
t2=time.time()

print('time',t2-t1)
print('UF:',R_uf)
print('subopt:',R_sub)
print('opt',R_aads)
print('mul_opt',R_mul)
print('div_opt',R_div)
print('mul_w',R_mul_w)
print('div_w',R_div_w)