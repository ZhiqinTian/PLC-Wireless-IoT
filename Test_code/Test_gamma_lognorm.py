import numpy as np
from numpy.core.numeric import Inf
from scipy.stats import lognorm, gamma, norm
from scipy.special import factorial, gammainc, erfcinv
from scipy.optimize import fsolve
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import time


m = 1
x1 = np.linspace(gamma.ppf(0.001, m, scale=1/m),
                gamma.ppf(0.999, m, scale=1/m), 100)
s = np.sqrt(np.log(1+1/m))
x2=np.linspace(lognorm.ppf(0.001, s, scale=1/np.exp(s**2/2)),lognorm.ppf(0.999, s, scale=1/np.exp(s**2/2)),100)

x=np.linspace(0,5,1000)

plt.plot(x, gamma.pdf(x, m, scale=1/m), '-',
         lw=1.5, alpha=0.6, label='gamma pdf')
plt.plot(x,lognorm.pdf(x, s, scale=1/np.sqrt(1+1/m)), '-',
         lw=1.5, alpha=0.6, label='lognorm pdf')

plt.legend(loc='best', frameon=False)
plt.xlim(0,5)
plt.show()
