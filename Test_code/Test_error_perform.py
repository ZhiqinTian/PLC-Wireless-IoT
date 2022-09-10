from Phy import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator

phy = PLC_Phy()

plt.figure(1)
plt.xlabel('sinr /dB')
plt.ylabel('ber')
plt.ylim(-5, 0)
plt.xlim(0, 31)
sinr = np.linspace(0, 31, 200)
color = {'BPSK': 'red', 'QPSK': 'green','QAM8':'saddlebrown', 'QAM16': 'blue','QAM32':'m', 'QAM64': 'black','QAM128': 'purple','QAM256': 'c'}
linestyle = {'1_2': '-', '2_3': '--', '3_4': ':'}
ax = plt.gca()
ax.xaxis.set_major_locator(MultipleLocator(1))
ax.xaxis.set_minor_locator(MultipleLocator(0.1))
ax.tick_params(axis='x',direction='in')
for m in MCS.get_Ms():
    for c in MCS.get_Cs():
        ber = np.array([phy.error_rate_model.ber_calc(
            MCS(m, c), pow(10, s/10)) for s in sinr])
        plt.plot(sinr, np.log10(ber),
                 color=color[m], linestyle=linestyle[c], label=m+'_'+c)
plt.legend()
plt.grid()
plt.show()
