from scipy.stats import lognorm, gamma, norm
from scipy.special import factorial, gammainc, erfcinv
from scipy.optimize import fsolve
import cvxpy as cv 
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np


# p_e_a=0.1
# N=200
#p_error = 1-np.power(1-p_e_a,1/N)
p_error=0.001
s = 1.2
m = 1.8


def pdf_hw_2(x):
    return gamma.pdf(x, m, scale=1 / m)

def pdf_hp_2(x):
    return lognorm.pdf(x, s, scale=1 / np.exp(s ** 2 / 2))

def cdf_hw_2(x):
    return gamma.cdf(x, m, scale=1 / m)

def cdf_hp_2(x):
    return lognorm.cdf(x, s, scale=1 / np.exp(s ** 2 / 2))+1e-8

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
    p_outage, _ = integrate.quad(f, 0, (2**R-1)/snr_p)
    return p_outage


# def R_max(p_outage, snr_w, m, snr_p, s):

#     R_w = np.log2(1+gamma.ppf(p_outage, m, scale=1/m)*snr_w)

#     R_p = np.log2(
#         1+np.exp(s*norm.ppf(p_outage)-s**2/2)*snr_p)

#     def sc_dive_R_max(log_R):
#         return np.log10(p_outage_sc_dive(2**log_R, snr_w, m, snr_p, s))-np.log10(p_outage)
#     R_sc = np.power(2, fsolve(sc_dive_R_max, [
#         np.log2(np.log2(1+max(snr_w, snr_p)))]))

#     def mrc_dive_R_max(log_R):
#         return np.log10(p_outage_mrc_dive(2**log_R, snr_w, m, snr_p, s))-np.log10(p_outage)
#     # R_mrc = np.power(2, fsolve(mrc_dive_R_max, [
#     #     np.log2(np.log2(1+snr_w+snr_p))]))
#     R_mrc = [0]
#     return R_w, R_p, R_sc[0], R_mrc[0]

def R_max(p_outage, snr_w, m, snr_p, s):
    
    R_w = np.log2(1+gamma.ppf(p_outage/2, m, scale=1/m)*snr_w)

    R_p = np.log2(
        1+np.exp(s*norm.ppf(p_outage/2)-s**2/2)*snr_p)

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

def R_mul_max(p_outage, snr_w, m, snr_p, s):
    
    
    R_w = np.log2(1+gamma.ppf(p_outage/2, m, scale=1/m)*snr_w)

    R_p = np.log2(1+np.exp(s*norm.ppf(p_outage/2)-s**2/2)*snr_p)
    
    def mul_R_max(log_R):
        f1=np.log(1-cdf_hw_2((2**(2**log_R[0])-1)/snr_w))+np.log(1-cdf_hp_2((2**(2**log_R[1])-1)/snr_p))-np.log(1-p_outage)
        #f1=(1-cdf_hw_2((2**(2**log_R[0])-1)/snr_w))*(1-cdf_hp_2((2**(2**log_R[1])-1)/snr_p))/(1-p_outage)
        f2=np.log(pdf_hw_2((2**(2**log_R[0])-1)/snr_w)*2**(2**log_R[0])/snr_w)-np.log(1-cdf_hw_2((2**(2**log_R[0])-1)/snr_w))+np.log(1-cdf_hp_2((2**(2**log_R[1])-1)/snr_p))-np.log(pdf_hp_2((2**(2**log_R[1])-1)/snr_p)*2**(2**log_R[1])/snr_p)
        # print(2**log_R)
        # print(f1,f2)
        return np.array([f1,f2])
    
    R_mul=np.power(2,fsolve(mul_R_max, [
        np.log2(R_w),np.log2(R_p)]))
    return np.sum(R_mul)

def R_mul_opt(p_outage, snr_w, snr_p):
    x = cv.Variable() # 定义变量x,定义变量y。两个都是标量
    y = cv.Variable()
    # Create two constraints.
    # 定义两个约束式
    constraints = [x>=0,y >= 0,np.log(1-cdf_hw_2((cv.exp(cv.multiply(np.log(2), x))-1)/snr_w))+np.log(1-cdf_hp_2((cv.exp(cv.multiply(np.log(2), y))-1)/snr_p))-np.log(1-p_outage)>=0]
    # 优化的目标函数
    obj = cv.Minimize(-x-y)
    # 把目标函数与约束传进Problem函数中
    prob = cv.Problem(obj, constraints)
    prob.solve()  # Returns the optimal value.
    return prob.value


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

#print(R_mul_max(0.001, 100, m, 100, s))


SNR = np.linspace(5,40, 200)
snr = np.power(10, SNR/10)
k_p = np.exp(s*norm.ppf(p_error/2)-s**2/2)
k_w = gamma.ppf(p_error/2, m, scale=1/m)

print('s:', s, 'm:', round(m, 2), 'k_p:', k_p, 'k_w:', k_w)
R_W = []
R_P = []
R_R = []

R_W_new = []
R_P_new = []
R_R_new = []

R_M=[]

R_SC = []
R_MRC = []
for snr_i in snr:
    #R_w, R_p, R_sc, R_mrc = R_max(p_error, snr_i, m, snr_i, s)
    R_w, R_p, R_sc, R_mrc = R_max(p_error, snr_i, m, snr_i, s)
    R_m=R_mul_max(p_error, snr_i, m, snr_i, s)
    
    R_W.append(R_w)
    R_P.append(R_p)
    R_M.append(R_m)
    # R_W_new.append(R_w_new)
    # R_P_new.append(R_p_new)
    
    R_SC.append(R_sc)
    R_MRC.append(R_mrc)
    R_R.append(R_w+R_p)
    #R_R_new.append(R_w_new+R_p_new)

plt.rc('font',family='Times New Roman')
plt.rc('text', usetex=True)

plt.figure(1,figsize=(7, 5), dpi=80)
plt.plot(SNR, R_SC, linestyle='--', label='Diversity',color='b')
#plt.plot(SNR, R_R, linestyle='-', label='Multiplexing',color='r')
plt.plot(SNR, R_R, linestyle='-', label='Multiplexing',color='r')
#plt.plot(SNR, R_M, linestyle=':', label='Multiplexing-opt',color='r')
# plt.plot(SNR, R_W, linestyle=':', label='w',color='r')
# plt.plot(SNR, R_P, linestyle=':', label='p',color='b')
plt.ylim(0, 10)
plt.xlim(5, 40)
plt.ylabel('Outage Rate (bps/Hz)',fontsize=17)
plt.xlabel('Average SNR (dB)',fontsize=17)
plt.grid(linestyle=':')
plt.rcParams.update({'font.size': 17})
plt.xticks(fontsize=15)
plt.yticks(fontsize=15) 
plt.legend()
plt.savefig("E:\\paper\\Paper\\fig\\Div_Mul.pdf",bbox_inches='tight',pad_inches = 0)
plt.figure(2,figsize=(5, 4), dpi=80)
r1_c=np.array([55.18,56.42,50.1,54.23,56.08,45.94,44.4,46.23,45.12,44.1])
r1=np.array([51.14,52.43,46.49,50.64,51.2,37.7,39.3,40.45,38.61,36.93])
r2_c=np.array([70.22,72.58,70.78,67.94,71.32,72.8,69.62,72.96,71.51,72.48])
r2=np.array([68.68,71,69.08,66.29,69.89,71.26,67.94,71.42,69.81,71.14])
n=np.linspace(1,10,10)
plt.plot(n, r1_c, linestyle='--', label='With compensation',color='r',marker='^',markerfacecolor='none',markersize=6)
plt.plot(n, r1, linestyle='-', label='Without compensation',color='b',marker='o',markerfacecolor='none',markersize=6)
plt.ylim(30, 80)
plt.xlim(0.7, 10.6)
plt.ylabel('Outage Rate (Mbps)',fontsize=20)
plt.xlabel('Samples',fontsize=20)
plt.grid(linestyle=':')
plt.rcParams.update({'font.size': 20})
plt.xticks(fontsize=20)
plt.yticks(fontsize=20) 
plt.legend()
#plt.savefig("D:\\code\\paper\\fig\\Po_case1.pdf",bbox_inches='tight',pad_inches = 0)
plt.savefig("E:\\paper\\Paper\\fig\\Po_case1.pdf",bbox_inches='tight',pad_inches = 0)

plt.figure(3,figsize=(5, 4), dpi=80)
plt.plot(n, r2_c, linestyle='--', label='With compensation',color='r',marker='^',markerfacecolor='none',markersize=6)
plt.plot(n, r2, linestyle='-', label='Without compensation',color='b',marker='o',markerfacecolor='none',markersize=6)
plt.ylim(55, 80)
plt.xlim(0.7, 10.6)
plt.ylabel('Outage Rate (Mbps)',fontsize=20)
plt.xlabel('Samples',fontsize=20)
plt.grid(linestyle=':')
plt.rcParams.update({'font.size': 20})
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend()
#plt.savefig("D:\\code\\paper\\fig\\Po_case2.pdf",bbox_inches='tight',pad_inches = 0)
#plt.savefig("E:\\paper\\Paper\\fig\\Po_case2.pdf",bbox_inches='tight',pad_inches = 0)
plt.figure(4,figsize=(7, 5), dpi=80)
SNR = np.linspace(5,40, 20)
snr = np.power(10, SNR/10)

SNR_p=np.array([10,30,40])
snr_p=np.power(10,SNR_p/10)

t=0
for r_p in snr_p:
    R_R = []
    R_M=[]
    for snr_i in snr:
        #R_w, R_p, R_sc, R_mrc = R_max(p_error, snr_i, m, snr_i, s)
        R_w, R_p, R_sc, R_mrc = R_max(p_error, snr_i, m, r_p, s)
        R_m=R_mul_max(p_error, snr_i,m,  r_p,s)
        R_M.append(R_m)
        R_R.append(R_w+R_p)
    #plt.plot(SNR,R_M)
    if t==0:
        plt.plot(SNR, R_R, linestyle='--', label='approximation',color='b')
        plt.plot(SNR, R_M, linestyle='-', linewidth=1, label='original',color='r')
        t+=1
    else:
        plt.plot(SNR, R_R, linestyle='--',color='b')
        plt.plot(SNR, R_M, linestyle='-', linewidth=1,color='r')

plt.annotate(r'$\bar{\gamma}_{p,i}$='+str(SNR_p[0])+' dB',xy=(25,1.8 ),xytext=(25, 1.8), fontsize=18)
plt.annotate(r'$\bar{\gamma}_{p,i}$='+str(SNR_p[1])+' dB',xy=(16,5),xytext=(16, 5), fontsize=18)
plt.annotate(r'$\bar{\gamma}_{p,i}$='+str(SNR_p[2])+' dB',xy=(12,7.5 ),xytext=(12, 7.5), fontsize=18)
plt.ylim(0, 12)
plt.xlim(5, 40)
plt.ylabel('Outage Rate (bps/Hz)',fontsize=17)
plt.xlabel(r'$\bar{\gamma}_{w,i}$ (dB)',fontsize=17)
plt.grid(linestyle=':')
plt.rcParams.update({'font.size': 17})
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.legend(fontsize=15)
plt.savefig("E:\\paper\\Paper\\fig\\approx_r_mul.pdf",bbox_inches='tight',pad_inches = 0)

# plt.figure(6)
# x=np.linspace(-1,1,100)
# plt.plot(x,np.log(1-cdf_hp_2(2**(2**x))))
plt.show()
