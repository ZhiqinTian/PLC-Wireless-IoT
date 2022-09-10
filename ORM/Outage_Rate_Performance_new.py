from cProfile import label
from scipy.stats import lognorm, gamma, norm,ncx2
from scipy.special import factorial, gammainc, erfcinv
from scipy.optimize import fsolve
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np


p_error=0.001
e_2=0.1

# def pdf_hw_2(x):
#     return gamma.pdf(x, m, scale=1 / m)
# def cdf_hw_2(x):
#     return gamma.cdf(x, m, scale=1 / m)
def pdf_hw_2(x,hf_2):
    return ncx2.pdf(x, 2, hf_2*2/e_2, scale=e_2/2)

# def pdfw_test(x,hf_2):
#     return ncx2.pdf(x*2/e_2, 2, hf_2*2/e_2)*2/e_2

def cdf_hw_2(x,hf_2):
    return ncx2.cdf(x, 2, hf_2*2/e_2, scale=e_2/2)

def ppf_hw_2(x,hf_2):
    return ncx2.ppf(x, 2, hf_2*2/e_2, scale=e_2/2)



def pdf_hp_2(x,s):
    return lognorm.pdf(x, s, scale=1 / np.exp(s ** 2 / 2))

def cdf_hp_2(x,s):
    return lognorm.cdf(x, s, scale=1 / np.exp(s ** 2 / 2))

def ppf_hp_2(x,s):
    return lognorm.ppf(x, s, scale=1 / np.exp(s ** 2 / 2))


x1=np.linspace(0.1,5,1000)

plt.figure(11)
plt.plot(x1,np.log(cdf_hw_2(x1,15)))
plt.figure(12)
x2=np.linspace(0.01,0.75,1000)
plt.plot(x2,np.log(cdf_hp_2(x2,1)))
plt.show()
exit(0)

def p_outage_wireless(R, snr_w, hf_2):
    return cdf_hw_2((2**R-1)/snr_w,hf_2)

def p_outage_plc(R, snr_p,s):
    return cdf_hp_2((2**R-1)/snr_p,s)

def p_outage_sc_dive(R, snr_w, hf_2, snr_p,s):
    #print(p_outage_wireless(R, snr_w, hf_2))
    p_outage = p_outage_wireless(R, snr_w, hf_2)*p_outage_plc(R, snr_p,s)
    return p_outage

def R_max(p_outage, snr_w, hf_2, snr_p,s):
    k1=ppf_hw_2(p_outage, hf_2)*pdf_hw_2(ppf_hw_2(p_outage, hf_2), hf_2)/p_outage
    k2=ppf_hp_2(p_outage, s)*pdf_hp_2(ppf_hp_2(p_outage, s),s)/p_outage
    
    kw=ppf_hw_2(k2/(k1+k2)*p_outage, hf_2)
    kp=ppf_hp_2(k1/(k1+k2)*p_outage, s)
    
    R_w = np.log2(1+kw*snr_w)
    R_p = np.log2(1+kp*snr_p)
    
    C_w=np.log2(1+snr_w)
    C_p=np.log2(1+snr_p)
    
    Ra=0
    Rb=C_w+C_p
    dealta_log_p=1
    while np.abs(dealta_log_p)>0.05:
        Rc=(Ra+Rb)/2
        dealta_log_p=np.log(p_outage_sc_dive(Rc, snr_w, hf_2, snr_p,s))-np.log(p_outage)
        if dealta_log_p<0:
            Ra=Rc
        else:
            Rb=Rc
    R_sc = Rc
    # def mrc_dive_R_max(log_R):
    #     return np.log10(p_outage_mrc_dive(2**log_R, snr_w, m, snr_p, s))-np.log10(p_outage)
    # R_mrc = np.power(2, fsolve(mrc_dive_R_max, [
    #     np.log2(np.log2(1+snr_w+snr_p))]))
    #R_mrc = [0]
    return R_w, R_p, R_sc

def R_mul_max(p_outage, snr_w, hf_2, snr_p,s):
    k1=ppf_hw_2(p_outage, hf_2)*pdf_hw_2(ppf_hw_2(p_outage, hf_2), hf_2)/p_outage
    k2=ppf_hp_2(p_outage,s)*pdf_hp_2(ppf_hp_2(p_outage,s),s)/p_outage
    rho1=k2/(k1+k2)
    rho2=1-rho1
    kw=ppf_hw_2(rho1*p_outage, hf_2)
    kp=ppf_hp_2(rho2*p_outage, s)
    #print(kw,kp)
    R_w = np.log2(1+kw*snr_w)
    R_p = np.log2(1+kp*snr_p)
    def mul_R_max1(log_R):
        f1=np.log(1-cdf_hw_2((2**(2**log_R[0])-1)/snr_w,hf_2))+np.log(1-cdf_hp_2((2**(2**log_R[1])-1)/snr_p))-np.log(1-p_outage)
        #f1=(1-cdf_hw_2((2**(2**log_R[0])-1)/snr_w,hf_2))*(1-cdf_hp_2((2**(2**log_R[1])-1)/snr_p))/(1-p_outage)
        #f2=np.log(pdf_hw_2((2**(2**log_R[0])-1)/snr_w,hf_2)*2**(2**log_R[0])/snr_w)-np.log(1-cdf_hw_2((2**(2**log_R[0])-1)/snr_w,hf_2))+np.log(1-cdf_hp_2((2**(2**log_R[1])-1)/snr_p))-np.log(pdf_hp_2((2**(2**log_R[1])-1)/snr_p)*2**(2**log_R[1])/snr_p)
        f2=np.log(pdf_hw_2((2**(2**log_R[0])-1)/snr_w,hf_2)*2**(2**log_R[0])/snr_w)-np.log(pdf_hp_2((2**(2**log_R[1])-1)/snr_p)*2**(2**log_R[1])/snr_p)
        # print(2**log_R)
        # print(f1,f2)
        return np.array([f1,f2])
    
    def mul_R_max2(snr_w,snr_p,hf_2,s):
        
        x=np.linspace(0,1,100)
        kw=ppf_hw_2(x*p_error, hf_2)
        kp=ppf_hp_2((1-x)*p_error, s)
        R_w = np.log2(1+kw*snr_w)
        R_p = np.log2(1+kp*snr_p)
        R_sum=R_w+R_p
        #index_max=np.argmax(R_sum)
        return np.max(R_sum)
    
    # while np.abs(dealta_log_p)>0.05:
    #     Rc=(Ra+Rb)/2
    #     dealta_log_p=np.log(p_outage_sc_dive(Rc, snr_w, hf_2, snr_p,s))-np.log(p_outage)
    #     if dealta_log_p<0:
    #         Ra=Rc
    #     else:
    #         Rb=Rc
    # R_w0=R_w
    # R_p0=R_p
    #if snr_w<=8 and snr_p>=19:
    return R_w,R_p,mul_R_max2(snr_w,snr_p,hf_2,s)
    # else:
    #     R_mul=np.power(2,fsolve(mul_R_max1, [
    #         np.log2(R_w0),np.log2(R_p0)]))
    # return R_w,R_p,np.sum(R_mul)


SNR = np.linspace(5,40, 20)
snr = np.power(10, SNR/10)

R_W = []
R_P = []
R_R = []

R_W_new = []
R_P_new = []
R_R_new = []

R_M=[]

R_SC = []
R_MRC = []
hf=snr/1000
for i in range(len(snr)):
    snr_i=snr[i]
    R_w, R_p, R_sc = R_max(p_error,snr_i, hf[i], snr_i,1)
    R_W.append(R_w)
    R_P.append(R_p)
    R_SC.append(R_sc)
    R_R.append(R_w+R_p)

plt.rc('font',family='Times New Roman')
plt.rc('text', usetex=True)

plt.figure(1,figsize=(7, 5), dpi=80)
plt.plot(SNR, R_SC, linestyle='--', label='Diversity',color='b')
plt.plot(SNR, np.log2(1+snr), linestyle='--', label='Shannon',color='k')
#plt.plot(SNR, R_R, linestyle='-', label='Multiplexing',color='r')
plt.plot(SNR, R_R, linestyle='-', label='Multiplexing',color='r')
#plt.plot(SNR, R_M, linestyle=':', label='Multiplexing-opt',color='r')
plt.plot(SNR, R_W, linestyle=':', label='w',color='r')
plt.plot(SNR, R_P, linestyle=':', label='p',color='b')
plt.ylim(0, 15)
plt.xlim(0, 40)
plt.ylabel('Outage Rate (bps/Hz)',fontsize=17)
plt.xlabel('Average SNR (dB)',fontsize=17)
plt.grid(linestyle=':')
plt.rcParams.update({'font.size': 17})
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.legend()

plt.figure(4,figsize=(7, 5), dpi=80)
SNR = np.linspace(0,40, 30)
snr = np.power(10, SNR/10)

SNR_p=np.array([0,20,30,40])
snr_p=np.power(10,SNR_p/10)

t=0
for r_p in snr_p:
    R_R = []
    R_M=[]
    for snr_i in snr:
        #R_w, R_p, R_sc, R_mrc = R_max(p_error, snr_i, m, snr_i, s)
        # R_w, R_p, R_sc, R_mrc = R_max(p_error, snr_i, 1, r_p, s)
        R_w, R_p,R_m=R_mul_max(p_error, snr_i,1,r_p,1)
        R_M.append(R_m)
        R_R.append(R_w+R_p)
    #plt.plot(SNR,R_M)
    if t==0:
        plt.plot(SNR, R_M, linestyle='-', linewidth=1, label='accurate',color='b')
        plt.plot(SNR, R_R, linestyle='--',linewidth=1, label='approximation',color='r')
        t+=1
    else:
        plt.plot(SNR, R_M, linestyle='-', linewidth=1,color='b')
        plt.plot(SNR, R_R, linestyle='--',linewidth=1,color='r')

plt.annotate(r'$\bar{\gamma}_{p,i}$='+str(SNR_p[0])+' dB',xy=(22,2.7 ),xytext=(22, 2.7), fontsize=18)
plt.annotate(r'$\bar{\gamma}_{p,i}$='+str(SNR_p[1])+' dB',xy=(16,5),xytext=(16, 5), fontsize=18)
plt.annotate(r'$\bar{\gamma}_{p,i}$='+str(SNR_p[2])+' dB',xy=(13.5,7 ),xytext=(13.5, 7), fontsize=18)
plt.annotate(r'$\bar{\gamma}_{p,i}$='+str(SNR_p[3])+' dB',xy=(9,10 ),xytext=(9, 10), fontsize=18)
plt.ylim(0,14)
plt.xlim(0, 40)
plt.ylabel('Outage Rate (bps/Hz)',fontsize=17)
plt.xlabel(r'$\bar{\gamma}_{w,i}$ (dB)',fontsize=17)
plt.grid(linestyle=':')
plt.rcParams.update({'font.size': 17})
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.legend(fontsize=15)
#plt.savefig("E:\\paper\\Paper\\fig\\approx_r_mul.pdf",bbox_inches='tight',pad_inches = 0)

# plt.figure(6)
# x=np.linspace(-1,1,100)
# plt.plot(x,np.log(1-cdf_hp_2(2**(2**x))))
plt.show()
