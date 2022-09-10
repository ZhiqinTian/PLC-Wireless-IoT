from Network import *
from Phy import *
from scipy.stats import lognorm
from scipy.special import factorial
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 1)
s = 1
mean, var, skew, kurt = lognorm.stats(s, moments='mvsk')
print(mean)
x = np.linspace(lognorm.ppf(0.01, s), lognorm.ppf(0.99, s), 100)
ax.plot(x, lognorm.pdf(x, s), 'r-', lw=5, alpha=0.6, label='lognorm pdf')
rv = lognorm(s)
ax.plot(x, rv.pdf(x), 'k-', lw=2, label='frozen pdf')
vals = lognorm.ppf([0.001, 0.5, 0.999], s)
np.allclose([0.001, 0.5, 0.999], lognorm.cdf(vals, s))
r = lognorm.rvs(s, loc=0, size=1000)
ax.hist(r, density=True, histtype='stepfilled', alpha=0.2)
ax.legend(loc='best', frameon=False)
h = lognorm.rvs(1,size=50)
t = np.linspace(1,20,50)
plt.show()
plt.plot(t,h)
plt.show()
R = np.linspace(1, 10, 100000)
for snr in range(4, 20):
    SNR = np.power(10, snr/10)
    #p_ray = np.log10(1-np.exp(-(np.power(2, R)-1)/SNR))
    #plt.plot(R, p_ray, label='SNR='+str(snr)+'dB')
    p_nakagami = np.log10(
        lognorm.cdf((np.power(2, R)-1)/SNR, s))
    plt.plot(R, p_nakagami, linestyle='-.')
plt.ylim(-6, 0.1)
plt.xlim(1, 10)
plt.xlabel('R b/s/Hz')
plt.ylabel('P_outage')
plt.legend()
plt.show()
# m = np.linspace(0,100,101)
# A = 3
# p = math.exp(-A)*np.power(A, m)/factorial(m)
# plt.plot(m,p)
# plt.show()
