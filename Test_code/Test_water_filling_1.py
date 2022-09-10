import numpy as np
from numpy.core.numeric import Inf
from scipy.stats import lognorm, gamma, norm
from scipy.special import factorial, gammainc, erfcinv
from scipy.optimize import fsolve
import scipy.integrate as integrate

from Network import *
from Phy import *
import matplotlib.pyplot as plt


s = 1.5
m = 1
p_error = 0.001


def p_pdf_hw_2(x):
    return gamma.pdf(x, m, scale=1 / m)


def p_pdf_hp_2(x):
    return lognorm.pdf(x, s, scale=1 / np.exp(s ** 2 / 2))


def p_pdf_hw_2_d(x):
    return ((m-1)/x-m)*gamma.pdf(x, m, scale=1 / m)


def p_pdf_hp_2_d(x):
    x1 = np.clip(x, a_min=0.001, a_max=None)
    return -(np.log(x1)/s**2+3)/x1*lognorm.pdf(x1, s, scale=1 / np.exp(s ** 2 / 2))


def p_cdf_hw_2(x):
    return gamma.cdf(x, m, scale=1 / m)


def p_cdf_hp_2(x):
    return lognorm.cdf(x, s, scale=1 / np.exp(s ** 2 / 2))


def p_outage_wireless(R, snr_w, m):
    return gamma.cdf((2 ** R - 1) / snr_w, m, scale=1 / m)


def p_outage_plc(R, snr_p, s):
    p_outage = lognorm.cdf((2 ** R - 1) / snr_p, s,
                           scale=1 / np.exp(s ** 2 / 2))
    return p_outage


def p_outage_sc_dive(R, snr_w, m, snr_p, s):
    p_outage = p_outage_wireless(R, snr_w, m) * p_outage_plc(R, snr_p, s)
    return p_outage


def p_outage_mrc_dive(R, snr_w, m, snr_p, s):
    def f(x):
        f = lognorm.pdf(x, s, scale=1 / np.exp(s ** 2 / 2)) * gamma.cdf(
            (2 ** R - 1 - x * snr_p) / snr_w, m, scale=1 / m
        )
        return f

    p_outage, _ = integrate.quad(f, 0, (2 ** R - 1) / snr_p)
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
    # R_mrc = np.power(2, fsolve(mrc_dive_R_max, [
    #     np.log2(np.log2(1+snr_w+snr_p))]))
    R_mrc = [0]
    return R_w, R_p, R_sc[0], R_mrc[0]


SNR = np.linspace(0, 50, 50)
snr = np.power(10, SNR/10)
R_W = []
R_P = []
R_R = []
R_SC = []
R_MRC = []
for snr_i in snr:
    R_w, R_p, R_sc, R_mrc = R_max(p_error, snr_i, m, snr_i, s)
    #R_W_A.append(w_R_average(p_error, snr_i, m))
    #print(R_w, np.log2(1+snr_i))
    R_W.append(R_w)
    R_P.append(R_p)
    R_SC.append(R_sc)
    R_MRC.append(R_mrc)
    R_R.append(R_w+R_p)
plt.figure(1)
snr_p = snr
snr_w = snr

plt.plot(SNR, R_W, linestyle='-.', label='wireless')
plt.plot(SNR, R_P, linestyle='--', label='plc')
plt.plot(SNR, R_SC, linestyle='--', label='sc')
#plt.plot(SNR, R_MRC, linestyle='--', label='mrc')
plt.plot(SNR, R_R, linestyle=':', label='reuse')
plt.ylim(0, 15)
plt.xlim(0, 50)
plt.ylabel('R b/s/Hz')
plt.xlabel('SNR')
plt.legend()

N = 100
sub_band = 20e6 / N
P_w_all = 100
P_p_all = 100
snr_w_ref = np.round(np.clip(np.random.normal(40, 10, size=N), 10, 50), 2)
snr_p_ref = np.round(np.clip(np.random.normal(40, 10, size=N), 10, 50), 2)
# plt.plot(np.linspace(1,N,N),snr_p_ref)

cnr_w = np.power(10, snr_w_ref / 10) * N / P_w_all
cnr_p = np.power(10, snr_p_ref / 10) * N / P_p_all


k_p = np.exp(s * norm.ppf(p_error) - s ** 2 / 2)
k_w = gamma.ppf(p_error, m, scale=1 / m)
factor_w = cnr_w * k_w
factor_p = cnr_p * k_p
P_w = np.clip((P_w_all + np.sum(1 / factor_w)) /
              N - 1 / factor_w, a_min=0.000001, a_max=None)
P_p = np.clip((P_p_all + np.sum(1 / factor_p)) / N -
              1 / factor_p, a_min=0.000001, a_max=None)
P_w = P_w * P_w_all / np.sum(P_w)
P_p = P_p * P_p_all / np.sum(P_p)
R_w = np.log2(1+factor_w*P_w)
R_p = np.log2(1+factor_p*P_p)
R_w_all = np.sum(R_w)
R_p_all = np.sum(R_p)
R_old = R_w_all+R_p_all
#print(R_w_all, R_p_all, R_w_all+R_p_all, sub_band*1e-6*(R_w_all+R_p_all))
#print('R_w:', R_w)
#print('R_p:', R_p)
p_sc = p_outage_sc_dive(R_w+R_p, P_w*cnr_w, m, P_p*cnr_p, s)
d = np.where(p_sc < p_error)[0]
r = np.where(p_sc >= p_error)[0]
print('d:', d)
print('r:', r)
Nd = len(d)
Nr = len(r)
factor_w_r = factor_w[r]
factor_p_r = factor_p[r]
cnr_w_d = cnr_w[d]
cnr_p_d = cnr_p[d]


def best_power_allocation(x):

    #x = np.clip(x, a_min=0.000001, a_max=None)
    
    P_w_d = x[0:Nd]
    P_p_d = x[Nd:2*Nd]
    
    snr_div = x[2*Nd:3*Nd]
    snr_w_d = P_w_d*cnr_w_d
    snr_p_d = P_p_d*cnr_p_d

    fw = p_pdf_hw_2(snr_div/snr_w_d)
    fp = p_pdf_hp_2(snr_div/snr_p_d)
    Fw = p_cdf_hw_2(snr_div/snr_w_d)
    Fp = p_cdf_hp_2(snr_div/snr_p_d)
    # print('Fw:',Fw,'Fp:',Fp,'fw:',fw,'fp:',fp)
    if (snr_div<=0).any():
        print('snr_div <=0')
    
    if (snr_w_d<=0).any():
        print('snr_w_d <=0')
    
    if (snr_p_d<=0).any():
        print('snr_p_d <=0')
    
    if (fw==0).any():
        print('fw 0')
        pw_error = list(filter(lambda x: x <= 0, P_w_d))
        print('pw error:',pw_error)
    
    if (fp==0).any():
        print('fp 0')
        pp_error = list(filter(lambda x: x <= 0, P_p_d))
        print('pp error:',pp_error)
    
    if (Fp==0).any():
        print('Fp 0')
        pp_error = list(filter(lambda x: x <= 0, P_p_d))
        print('pp error:',pp_error)
    
    if (Fw==0).any():
        print('Fw 0')
        pw_error = list(filter(lambda x: x <= 0, P_w_d))
        print('pw error:',pw_error)
    delta_w = np.sum(P_w_d)-P_w_all-np.sum(1/factor_w_r)
    delta_p = np.sum(P_p_d)-P_p_all-np.sum(1/factor_p_r)
    yw = P_w_d*(1+1/snr_div)*(1+snr_w_d * fp*Fw/(snr_p_d*fw*Fp))+delta_w/Nr

    yp = P_p_d*(1+1/snr_div)*(1+snr_p_d * fw*Fp/(snr_w_d*fp*Fw))+delta_p/Nr

    # yw = P_w_d * (1+snr_div)*(snr_w_d * fp*Fw+snr_p_d*fw*Fp) + \
    #     snr_div*snr_w_d*fw*Fp*delta_w / Nr
    # yp = P_p_d * (1+snr_div)*(snr_p_d * fw*Fp+snr_w_d*fp*Fw) + \
    #     snr_div*snr_p_d*fp*Fw*delta_p / Nr

    log_pe = np.log(p_error)-np.log(1e-12+Fp)-np.log(1e-12+Fw)

    result = np.concatenate((yw, yp, log_pe), axis=0)
    return result


x_w = np.ones(Nd)*P_w_all/N
x_p = np.ones(Nd)*P_p_all/N
# x_w = P_w[d]
# x_p = P_p[d]
x0 = np.concatenate(
    (x_w, x_p, np.power(2, R_w[d]+R_p[d]+0.00001)-1), axis=0)
result = np.round(fsolve(best_power_allocation, x0), 2)
#result = x0
#print('x0', np.round(x0, 2))

P_w_d = result[0:Nd]
P_p_d = result[Nd:2*Nd]
snr_div = result[2*Nd:3*Nd]

P_w_r = (P_w_all-np.sum(P_w_d)+np.sum(1/factor_w_r))/Nr - 1/factor_w_r
P_p_r = (P_p_all-np.sum(P_p_d)+np.sum(1/factor_p_r))/Nr - 1/factor_p_r
R_d = np.log2(1+snr_div)
R_r = np.log2(1+factor_w_r*P_w_r) + np.log2(1+factor_p_r*P_p_r)
R_new = np.sum(R_d)+np.sum(R_r)
#print('P_w_d', P_w_d, 'P_p_d', P_p_d, 'C_d', R_d)
print('old:',R_old,'new:',R_new)
plt.figure(2)
plt.bar(np.linspace(1, N, N), snr_p_ref)
plt.xlabel('subcarriers')
plt.ylabel('snr_p_ref')
plt.figure(3)
plt.bar(np.linspace(1, N, N), p_sc)
plt.xlabel('subcarriers')
plt.ylabel('p_sc')
plt.figure(4)
plt.bar(np.linspace(1, N, N), P_w)
plt.xlabel('subcarriers')
plt.ylabel('P_w')
plt.figure(5)
plt.bar(np.linspace(1, N, N), snr_w_ref)
plt.xlabel('subcarriers')
plt.ylabel('snr_w_ref')
plt.figure(6)
plt.bar(np.linspace(1, N, N), 10*np.log10(P_w*cnr_w))
plt.xlabel('subcarriers')
plt.ylabel('snr_w')
plt.figure(7)
plt.bar(np.linspace(1, N, N), 10*np.log10(P_p*cnr_p))
plt.xlabel('subcarriers')
plt.ylabel('snr_p')
plt.figure(8)
plt.bar(np.linspace(1, N, N), P_p)
plt.xlabel('subcarriers')
plt.ylabel('P_p')
# plt.ylim(0, 1.5)
#plt.xlim(0, 100)
plt.show()
