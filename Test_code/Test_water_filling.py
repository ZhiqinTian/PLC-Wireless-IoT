import numpy as np
from numpy.core.numeric import Inf
from scipy.stats import lognorm, gamma, norm
from scipy.special import factorial, gammainc, erfcinv
from scipy.optimize import fsolve
import scipy.integrate as integrate

from Network import *
from Phy import *
import matplotlib.pyplot as plt

# nodes = Node.creat_nodes(10)
# topo = Topology()
# pos = [
#     [0, 0],
#     [5, 0],
#     [8, 0],
#     [120, 0],
#     [5, 5],
#     [8, 3],
#     [5, -5],
#     [15, 0],
#     [15, 4],
#     [15, -3],
# ]
# topo.set_nodes(nodes, pos)

# topo.add_line(nodes[0], nodes[1])
# topo.add_line(nodes[1], nodes[2])
# topo.add_line(nodes[2], nodes[7])
# topo.add_line(nodes[7], nodes[3])
# topo.add_line(nodes[1], nodes[4])
# topo.add_line(nodes[2], nodes[5])
# topo.add_line(nodes[1], nodes[6])
# topo.add_line(nodes[7], nodes[8])
# topo.add_line(nodes[7], nodes[9])
# topo.show()

# plc_sm = PLC_Phy.Spectrum_Model(2e6, 30e6, 110)
# plc_channel = PLC_Channel(plc_sm)
# plc_channel.set_topo(topo)
# plc_channel.nodes_transfers_calc(nodes)

# rf_sm = RF_Phy.Spectrum_Model(2.4e9, 2e6, 30e6, 110)
# rf_channel = RF_Channel(rf_sm)
# rf_channel.set_topo(topo)

# plc_phy = PLC_Phy()
# rf_phy = RF_Phy()

# f = np.array(plc_sm.get_subcarriers())
# transfer1 = plc_channel.get_transfer(nodes[0], nodes[3])
# transfer2 = rf_channel.pathloss_calc(nodes[0], nodes[3])

# plc_tx_psd = Link_Adap.ave_psd_alloc(plc_sm, 100)
# # plc_noise_psd = np.power(
# #     10, (-145 + 53.23 * np.power(plc_sm.subcarriers * 1e-6, -0.337)) / 10
# # )
# plc_noise_psd = np.power(
#     10, (-140 + 53.23 * np.power(plc_sm.subcarriers * 1e-6, -0.337)) / 10
# )
# plc_rx_psd = plc_tx_psd * np.square(np.power(10, transfer1 / 10))
# plc_snr = plc_phy.sinr_model.snr_calc(plc_rx_psd, plc_noise_psd)

# rf_tx_psd = Link_Adap.ave_psd_alloc(rf_sm, 100)
# rf_noise_psd = np.array([4e-18] * rf_sm.num_subcarriers)
# rf_rx_psd = rf_tx_psd * np.square(np.power(10, transfer2 / 10))
# rf_snr = rf_phy.sinr_model.snr_calc(rf_rx_psd, rf_noise_psd)
# plt.figure(1)
# plt.plot(f, transfer1, label="PLC")
# plt.plot(f, transfer2, label="RF")
# plt.legend()
# plt.figure(2)
# plt.plot(f, 10 * np.log10(plc_snr), label="PLC")
# plt.plot(f, 10 * np.log10(rf_snr), label="RF")
# plt.legend()
# plt.show()

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
# mu_mrc = np.log(2*snr)-np.log(np.sqrt(1+(np.exp(s**2-1)+s_w)/4))
# sigma_mrc = np.sqrt(np.log(1+(np.exp(s**2-1)+s_w)/4))
#R_ap = np.log2(1+np.exp(norm.ppf(p_error)*sigma_mrc+mu_mrc))
#plt.plot(SNR, R_ap, linestyle=':', label='mrc_ap1')

#plt.plot(SNR, np.log2(1+k_p*snr_p+k_w*snr_w), linestyle=' ',marker = "x", label='mrc_ap2')
#plt.plot(SNR, np.log2(1+snr_p+snr_w), linestyle=':',marker = "+", label='mrc_ap3')
#plt.plot(SNR, np.array(R_MRC)-np.log2(1+k_p*snr_p+k_w*snr_w), linestyle=':', label='di')
#plt.plot(SNR, np.log2(1+snr_w), linestyle=':', label='ap3')
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

N = 20
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
              N - 1 / factor_w, a_min=0.0001, a_max=None)
P_p = np.clip((P_p_all + np.sum(1 / factor_p)) / N -
              1 / factor_p, a_min=0.0001, a_max=None)
P_w = P_w * P_w_all / np.sum(P_w)
P_p = P_p * P_p_all / np.sum(P_p)
R_w = np.log2(1+factor_w*P_w)
R_p = np.log2(1+factor_p*P_p)
R_w_all = np.sum(R_w)
R_p_all = np.sum(R_p)
print(R_w_all, R_p_all, R_w_all+R_p_all, sub_band*1e-6*(R_w_all+R_p_all))
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

    P_w_d = x[0:Nd]
    P_p_d = x[Nd:2*Nd]
    snr_div = x[2*Nd:3*Nd]
    snr_w_d = P_w_d*cnr_w_d
    snr_p_d = P_p_d*cnr_p_d

    fw = p_pdf_hw_2(snr_div/snr_w_d)
    fp = p_pdf_hp_2(snr_div/snr_p_d)
    Fw = p_cdf_hw_2(snr_div/snr_w_d)
    Fp = p_cdf_hp_2(snr_div/snr_p_d)

    delta_w = np.sum(P_w_d)-P_w_all-np.sum(1/factor_w_r)
    delta_p = np.sum(P_p_d)-P_p_all-np.sum(1/factor_p_r)
    # yw = snr_div*snr_w_d*fw*Fp*(P_w_all - np.sum(P_w_d)+np.sum(1/factor_w_r)) / \
    #     Nr-P_w_d * (1+snr_div)*(snr_w_d * fp*Fw+snr_p_d*fw*Fp)
    # yp = snr_div*snr_p_d*fp*Fw*(P_p_all - np.sum(P_p_d)+np.sum(1/factor_p_r)) / \
    #     Nr-P_p_d * (1+snr_div)*(snr_p_d * fw*Fp+snr_w_d*fp*Fw)

    yw = P_w_d * (1+snr_div)*(snr_w_d * fp*Fw+snr_p_d*fw*Fp) + \
        snr_div*snr_w_d*fw*Fp*delta_w / Nr
    yp = P_p_d * (1+snr_div)*(snr_p_d * fw*Fp+snr_w_d*fp*Fw) + \
        snr_div*snr_p_d*fp*Fw*delta_p / Nr

    log_pe = np.log(p_error)-np.log(1e-10+Fp)-np.log(1e-10+Fw)

    result = np.concatenate((yw, yp, log_pe), axis=0)
    return result


def jac(x):

    P_w_d = x[0:Nd]
    P_p_d = x[Nd:2*Nd]
    snr_div = x[2*Nd:3*Nd]
    snr_w_d = P_w_d*cnr_w_d
    snr_p_d = P_p_d*cnr_p_d

    fw = p_pdf_hw_2(snr_div/snr_w_d)
    fp = p_pdf_hp_2(snr_div/snr_p_d)
    Fw = p_cdf_hw_2(snr_div/snr_w_d)
    Fp = p_cdf_hp_2(snr_div/snr_p_d)

    fw_d = p_pdf_hw_2_d(snr_div/snr_w_d)
    fp_d = p_pdf_hp_2_d(snr_div/snr_p_d)

    delta_w = np.sum(P_w_d)-P_w_all-np.sum(1/factor_w_r)
    delta_p = np.sum(P_p_d)-P_p_all-np.sum(1/factor_p_r)

    yw_w_diff_diag = (1+snr_div)*(2*snr_w_d*fp*Fw+snr_p_d*fw*Fp-snr_div*fp*fw-snr_div*(
        snr_p_d/snr_w_d)*fw_d*Fp)+snr_div*snr_p_d*Fp*(fw-snr_div/(P_w_d*snr_w_d)*fw_d*delta_w)/Nr
    yw_w_diff_non_diag = (snr_div*snr_p_d*fw*Fp)/Nr

    yp_p_diff_diag = (1+snr_div)*(2*snr_p_d*fw*Fp+snr_w_d*fp*Fw-snr_div*fw*fp-snr_div*(
        snr_w_d/snr_p_d)*fp_d*Fw)+snr_div*snr_w_d*Fw*(fp-snr_div/(P_p_d*snr_p_d)*fp_d*delta_p)/Nr
    yp_p_diff_non_diag = (snr_div*snr_w_d*fp*Fw)/Nr

    yw_p_diff_diag = (1+snr_div)*(P_w_d*cnr_p_d*fw*Fp-snr_div*P_w_d*snr_w_d/(P_p_d*snr_p_d) *
                                  fp_d*Fw-snr_div*P_w_d/P_p_d*fw*fp)+snr_div*(cnr_p_d*fw*Fp-snr_div/P_p_d*fw*fp)*delta_w/Nr
    yp_w_diff_diag = (1+snr_div)*(P_p_d*cnr_w_d*fp*Fw-snr_div*P_p_d*snr_p_d/(P_w_d*snr_w_d) *
                                  fw_d*Fp-snr_div*P_p_d/P_w_d*fp*fw)+snr_div*(cnr_w_d*fp*Fw-snr_div/P_w_d*fp*fw)*delta_p/Nr

    yw_div_diff_diag = P_w_d*(snr_w_d*fp*Fw+snr_p_d*fw*Fp+(1+snr_div)*(2*fp*fw+snr_w_d/snr_p_d*fp_d*Fw +
                              snr_p_d/snr_w_d*fw_d*Fp))+(snr_p_d*fw*Fp+snr_div*fw*fp+snr_div*snr_p_d/snr_w_d*fp_d*Fw)*delta_w/Nr
    yp_div_diff_diag = P_p_d*(snr_p_d*fw*Fp+snr_w_d*fp*Fw+(1+snr_div)*(2*fw*fp+snr_p_d/snr_w_d*fw_d*Fp +
                              snr_w_d/snr_p_d*fp_d*Fw))+(snr_w_d*fp*Fw+snr_div*fp*fw+snr_div*snr_w_d/snr_p_d*fw_d*Fp)*delta_p/Nr

    y_sc_w_diff_diag = -fw*snr_div/((Fw+1e-12)*P_w_d*snr_w_d)
    y_sc_p_diff_diag = -fp*snr_div/((Fp+1e-12)*P_p_d*snr_p_d)
    y_sc_div_diff_diag = fw/((Fw+1e-12)*snr_w_d) + fp/((Fp+1e-12)*snr_p_d)

    J = np.zeros((3*Nd, 3*Nd))
    for i in range(Nd):
        J[i, 0:Nd] = yw_w_diff_non_diag[i]
        J[i, i] = yw_w_diff_diag[i]
        J[i, i+Nd] = yw_p_diff_diag[i]
        J[i, i+2*Nd] = yw_div_diff_diag[i]

        J[i+Nd, i] = yp_w_diff_diag[i]
        J[i+Nd, 0:Nd] = yp_p_diff_non_diag[i]
        J[i+Nd, i+Nd] = yp_p_diff_diag[i]
        J[i+Nd, i+2*Nd] = yp_div_diff_diag[i]

        J[i+2*Nd, i] = y_sc_w_diff_diag[i]
        J[i+2*Nd, i+Nd] = y_sc_p_diff_diag[i]
        J[i+2*Nd, i+2*Nd] = y_sc_div_diff_diag[i]
    return J


# x_w = np.ones(Nd)*P_w_all/N
# x_p = np.ones(Nd)*P_p_all/N
x_w = P_w[d]
x_p = P_p[d]
x0 = np.concatenate(
    (x_w, x_p, np.power(2, R_w[d]+R_p[d]+0.001)-1), axis=0)
result = np.round(fsolve(best_power_allocation, x0, fprime=jac), 2)
#result = x0
print('x0', np.round(x0, 2))
print(result)

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
