import numpy as np
from numpy.core.numeric import Inf
from scipy.stats import lognorm, gamma, norm
from scipy.special import factorial, gammainc, erfcinv
from scipy.optimize import fsolve
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import time


def p_outage_wireless(R, snr_w, m):
    return gamma.cdf((2**R-1)/snr_w, m, scale=1/m)


def p_outage_plc(R, snr_p, s):
    p_outage = lognorm.cdf((2**R-1)/snr_p, s, scale=1/np.exp(s**2/2))
    return p_outage


def p_outage_sc_dive(R, snr_w, m, snr_p, s):
    p_outage = p_outage_wireless(R, snr_w, m)*p_outage_plc(R, snr_p, s)
    return p_outage


def p_outage_mrc_dive(R, snr_w, m, snr_p, s):

    def f(x):
        f = lognorm.pdf(x, s, scale=1/np.exp(s**2/2)) * \
            gamma.cdf((2**R-1-x*snr_p)/snr_w, m, scale=1/m)
        return f
    p_outage, _ = integrate.quad(f, 0, (2**R-1)/snr_p,epsabs=10e-1,epsrel=10e-1)
    return p_outage


def R_max(p_outage, snr_w, m, snr_p, s):

    R_w = np.log2(1+gamma.ppf(p_outage, m, scale=1/m)*snr_w)

    R_p = np.log2(
        1+np.exp(s*norm.ppf(p_outage)-s**2/2)*snr_p)

    def sc_dive_R_max(log_R):
        return np.log10(p_outage_sc_dive(2**log_R, snr_w, m, snr_p, s))-np.log10(p_outage)
    R_sc = np.power(2, fsolve(sc_dive_R_max, [
        np.log2(np.log2(1+max(snr_w, snr_p)))]))

    def mrc_dive_R_max(log_R):
        return np.log10(p_outage_mrc_dive(2**log_R, snr_w, m, snr_p, s))-np.log10(p_outage)
    R_mrc = np.power(2, fsolve(mrc_dive_R_max, [
        np.log2(np.log2(1+snr_w+snr_p))]))
    #R_mrc = [0]
    return R_w, R_p, R_sc[0], R_mrc[0]


def R(p_error, snr, n):
    C = np.log2(1+snr)
    V = snr*(snr+2)/((snr+1)**2)*(np.log2(np.exp(1))**2)
    R = C - np.sqrt(V/n)*erfcinv(2*p_error)
    return R


def w_R_average(p_error, snr_w, m):
    n = 1000

    def R_average(x):
        return gamma.pdf(x, m, scale=1/m)*R(p_error, snr_w*x, n)

    R_a, _ = integrate.quad(R_average, 0, 20)
    return R_a


SNR = np.linspace(0, 30, 50)
snr = np.power(10, SNR/10)
s_p = 1.5
print(np.e**(s_p**2)-1)
s_w = 0.9
m = 1/(np.power(2*s_w+1, s_w)-1)
p_error = 0.01
s = s_p
m = 1
# M =np.linspace(1,5,10000)
# S =np.linspace(0.001,2,1000)
# K1 = []
# K2 = []
# for m in M:
#     k_w = gamma.ppf(p_error, m, scale=1/m)
#     K1.append(k_w)
# for s in S:
#     #print(s,s*norm.ppf(p_error))
#     k_p = np.exp(s*norm.ppf(p_error)-s**2/2)
#     K2.append(k_p)
# plt.plot(np.log2(3+M),K1)
# #plt.plot(np.exp(-S),K2)
# plt.show()
# exit(0)
k_p = np.exp(s*norm.ppf(p_error)-s**2/2)
k_w = gamma.ppf(p_error, m, scale=1/m)
print('s:',s,'m:',round(m,2),'k_p:',k_p,'k_w:',k_w)
R_W = []
R_P = []
R_R = []
R_SC = []
R_MRC = []
#R_W_A = []

for snr_i in snr:
    t1 = time.time()
    R_w, R_p, R_sc, R_mrc = R_max(p_error, snr_i, m, 30, s)
    t2 = time.time()
    print('t:',t2-t1)
    #R_W_A.append(w_R_average(p_error, snr_i, m))
    #print(R_w, np.log2(1+snr_i))
    R_W.append(R_w)
    R_P.append(R_p)
    R_SC.append(R_sc)
    R_MRC.append(R_mrc)
    R_R.append(R_w+R_p)
R_W=np.array(R_W)
R_P=np.array(R_P)
R_SC=np.array(R_SC)
R_MRC=np.array(R_MRC)
R_R=np.array(R_R)
plt.figure(1)
plt.plot(snr, R_W, linestyle='-.', label='wireless')
plt.plot(snr, R_P, linestyle='--', label='plc')
plt.plot(snr, R_SC-R_P, linestyle='--', label='sc')
plt.plot(snr, R_MRC-R_P, linestyle='--', label='mrc')
#plt.plot(snr, 2**R_R-1, linestyle=':', label='reuse')
#plt.ylim(0, 15)
#plt.xlim(0, 50)
plt.ylabel('R b/s/Hz')
plt.xlabel('SNR')
plt.legend()
plt.show()
