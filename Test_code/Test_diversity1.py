from Network import *
from Phy import *
from scipy.stats import nakagami, gamma, lognorm
from scipy.special import factorial, erfcinv, erfc
import matplotlib.pyplot as plt

plt.figure(1)
R = np.linspace(0, 10, 100000)
s_p = 0.1
s_w = 0.1
m = 1/(np.power(2*s_w+1, s_w)-1)
snr_p = 30
snr_w = 30
SNR_P = np.power(10, snr_p/10)
SNR_W = np.power(10, snr_w/10)
#p_ray = np.log10(1-np.exp(-(np.power(2, R)-1)/SNR))
#plt.plot(R, p_ray, label='SNR='+str(snr)+'dB')
for s_p in np.linspace(0.1, 1, 10):
    p_lognorm = np.log10(lognorm.cdf((np.power(2, R)-1)/SNR_P, 2.6*s_p))
    plt.plot(np.log(R), p_lognorm, linestyle='--', label='s='+str(round(2.6*s_p, 1)))
plt.ylim(-10, 0.1)
plt.xlim(0, 10)
plt.xlabel('R b/s/Hz')
plt.ylabel('P_outage')
plt.legend()
plt.figure(5)
for s_w in np.linspace(0.1, 1, 10):
    m = 1/(np.power(2*s_w+1, s_w)-1)
    p_nakagami = np.log10(
        gamma.cdf((np.power(2, R)-1)/SNR_W, m, loc=0, scale=1/m))
    # p_nakagami = np.log10(
    #     gamma.cdf(R/SNR_W, m, loc=0, scale=1/m))
    plt.plot(R, p_nakagami, linestyle='-.', label='m='+str(round(m, 1)))
# p_diversity = p_lognorm+p_nakagami
# plt.plot(R, p_diversity, linestyle=':',label='diversity')
plt.ylim(-10, 0.1)
plt.xlim(0, 10)
plt.xlabel('R b/s/Hz')
plt.ylabel('P_outage')
plt.legend()

# plt.figure(2)
# SNR = np.linspace(-5, 20, 1000)
# p = np.power(10, SNR/10)
# R = np.linspace(0.5, 10, 20)
# n = 1000
# C = np.log2(1+p)
# v = p*(p+2)/((p+1)**2)*(np.log2(np.exp(1))**2)
# for r in R:
#     pe = 1/2*erfc((C - r)/np.sqrt(v/n))
#     plt.plot(SNR, np.log10(pe), label='R='+str(r))
# plt.ylim(-10, 0.1)
# plt.xlim(-5, 20)
# plt.xlabel('SNR dB')
# plt.ylabel('Pe')
# plt.legend()

# plt.figure(3)
# SNR = np.linspace(-5, 30, 18)
# n = 1000
# R = np.linspace(0.1, 12, 1000)
# for snr in SNR:
#     p = np.power(10, snr/10)
#     C = np.log2(1+p)
#     v = p*(p+2)/((p+1)**2)*(np.log2(np.exp(1))**2)
#     pe = 1/2*erfc((C - R)/np.sqrt(v/n))
#     plt.plot(R, np.log10(pe), label='SNR='+str(int(snr)))
#     plt.plot(np.ones(len(R))*C, np.log10(pe), linestyle='-.')
# plt.ylim(-10, 0.1)
# plt.xlim(0.1, 12)
# plt.xlabel('R b/s/Hz')
# plt.ylabel('Pe')
# plt.legend()
plt.show()
