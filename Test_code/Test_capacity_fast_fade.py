import numpy as np
from numpy.core.numeric import Inf
from scipy.stats import lognorm, gamma, norm
from scipy.special import factorial, gammainc, erfcinv
from scipy.optimize import fsolve
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import time


def C_w(snr_w, m):
    def f(x):
        f = gamma.pdf(x, m, scale=1 / m) * np.log2(1 + x * snr_w)
        return f

    C, _ = integrate.quad(f, 0, 20)
    return C


def C_p(snr_p, s):
    def f(x):
        f = lognorm.pdf(x, s, scale=1 / np.exp(s ** 2 / 2)) * np.log2(1 + x * snr_p)
        return f

    C, _ = integrate.quad(f, 0, 20)
    return C


def C_mrc(snr_w, m, snr_p, s):
    def f(x, y):
        f = (
            gamma.pdf(x, m, scale=1 / m)
            * lognorm.pdf(y, s, scale=1 / np.exp(s ** 2 / 2))
            * np.log2(1 + x * snr_w + y * snr_p)
        )
        return f

    C, _ = integrate.dblquad(f, 0, 20, 0, 20)
    return C


C1 = C_p(2, 1)
C2 = C_w(2, 1)
print(C1, C2, C_mrc(2, 1, 2, 1), C1 + C2, np.log2(1 + 4), 2 * np.log2(1 + 2))

SNR = np.linspace(1, 50, 50)
snr = np.power(10, SNR / 10)
s_p = 0.8
s_w = 0.8
p_error = 1e-4
m = 1 / (np.power(2 * s_w + 1, s_w) - 1)
s = (1.6 * s_p) ** (1 / s_p)


# plt.plot(SNR, R_W_A, linestyle='--', label='w_R_A')
# plt.plot(SNR, np.log2(1+0.017*snr), linestyle=':', label='ap1')
# plt.plot(SNR, np.log2(1+0.1*snr), linestyle=':', label='ap3')
# plt.plot(SNR, np.log2(1+snr), linestyle='--', label='shannon C')
# plt.plot(SNR, R_W, linestyle='-.', label='wireless')
# plt.plot(SNR, R_P, linestyle='--', label='plc')
# plt.plot(SNR, R_SC, linestyle='--', label='sc')
# plt.plot(SNR, R_MRC, linestyle='--', label='mrc')
# plt.plot(SNR, R_R, linestyle=':', label='reuse')
# plt.ylim(0, 20)
# plt.xlim(0, 50)
# plt.ylabel('R b/s/Hz')
# plt.xlabel('SNR')
# plt.legend()
# plt.show()
