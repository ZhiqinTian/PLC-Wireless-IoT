from scipy.stats import lognorm, gamma, norm
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import numpy as np

p_error = 0.001
s = 1.5
m = 1.2

def pdf_hw_2(x):
    return gamma.pdf(x, m, scale=1 / m)

def pdf_hp_2(x):
    return lognorm.pdf(x, s, scale=1 / np.exp(s ** 2 / 2))

def cdf_hw_2(x):
    return gamma.cdf(x, m, scale=1 / m)

def cdf_hp_2(x):
    return lognorm.cdf(x, s, scale=1 / np.exp(s ** 2 / 2))

def func_w(x, a, b):
    return 1-np.exp(-a*x**b)

def func_p(x, a, b):
    return np.exp(-a/x**b)

def p_error_div(R, snr_w, snr_p):
    p_error = cdf_hw_2((2**R-1)/snr_w)*cdf_hp_2((2**R-1)/snr_p)
    return p_error


def outage_R_div(snr_w, snr_p, p_error):
    def otage_R_div_eq(log_R):
        return np.log10(p_error_div(2**log_R, snr_w, snr_p))-np.log10(p_error)
    outage_R = np.power(2, fsolve(otage_R_div_eq, [
        np.log2(np.log2(1+max(snr_w, snr_p)))]))
    return outage_R[0]


def water_filling(omiga, p_all):
    N = len(omiga)
    flag = np.ones(N)
    p = np.ones(N)
    while True:
        flag[np.where(p < 0)] = 0
        miu = np.sum(flag)/(p_all+np.sum((1/omiga)*flag))
        p = (1/miu-1/omiga)*flag
        if (p >= 0).all():
            break
    return np.abs(p)


def water_filling_div_mul(omiga_d, omiga_m, p_all):
    Nm = len(omiga_m)
    flag_m = np.ones(Nm)
    p_m = np.ones(Nm)
    while True:
        flag_m[np.where(p_m < 0)] = 0
        miu = (np.sum(flag_m)+np.sum(omiga_d))/(p_all+np.sum((1/omiga_m)*flag_m))
        p_d = omiga_d/miu
        p_m = (1/miu-1/omiga_m)*flag_m
        if (p_m >= 0).all():
            break
    return p_d, np.abs(p_m)

def duality_water_filling(div,mul,k1,k2,k3,k_w_m,k_p_m,p_w_all,p_p_all):
    Nd=len(div)
    Nm=len(mul)
    p_w_new=np.ones(Nd+Nm)
    p_p_new=np.ones(Nd+Nm)
    v_d=0.01*np.ones(Nd)
    k_v=0.01
    while True:
        #print('v:',v_d)
        p_w_new[div],p_w_new[mul] = water_filling_div_mul(k1+v_d,k_w_m,p_w_all)
        p_p_new[div],p_p_new[mul] = water_filling_div_mul(k2+v_d,k_p_m,p_p_all)
        p_w_new=np.clip(p_w_new, a_min=0.001, a_max=None)
        p_p_new=np.clip(p_p_new, a_min=0.001, a_max=None)
        #print('contr:',k1*np.log(p_w_new[div])+k2*np.log(p_p_new[div])+k3)
        p_w_new=np.clip(p_w_new, a_min=0.001, a_max=None)
        p_p_new=np.clip(p_p_new, a_min=0.001, a_max=None)
        v_d_new = np.clip(v_d-k_v*(k1*np.log(p_w_new[div])+k2*np.log(p_p_new[div])+k3),a_min=0,a_max=None)
        if np.linalg.norm(v_d_new-v_d)<0.005:
            return p_w_new,p_p_new
        else:
            v_d = v_d_new

def get_approx_parameter(fi, p_w_d, p_p_d, cnr_w_d, cnr_p_d):
    x_w=fi/(p_w_d*cnr_w_d)
    x_p=fi/(p_p_d*cnr_p_d)
    fw=pdf_hw_2(x_w)
    Fw=np.clip(cdf_hw_2(x_w),a_min=1e-8,a_max=None)
    fp=pdf_hp_2(x_p)
    Fp=np.clip(cdf_hp_2(x_p),a_min=1e-8,a_max=None)
    alpha1=fw*x_w/Fw
    alpha2=fp*x_p/Fp
    # print('alpha1',alpha1)
    # print('alpha2',alpha2)
    beta1 = np.log(Fw)-alpha1*np.log(x_w)
    beta2 = np.log(Fp)-alpha2*np.log(x_p)
    omiga1=alpha1/(alpha1+alpha2)
    omiga2=alpha2/(alpha1+alpha2)
    omiga3=omiga1*np.log(cnr_w_d)+omiga2*np.log(cnr_p_d)-(beta1+beta2)/(alpha1+alpha2)+np.log(p_error)/(alpha1+alpha2)
    k0=fi/(fi+1)
    k1=omiga1*k0
    k2=omiga2*k0
    k3=omiga3*k0+np.log(1+fi)-k0*np.log(fi)
    #print('fi',fi)
    #print('fi_ap:',np.exp(k1*np.log(p_w_d)+k2*np.log(p_p_d)+k3)-1)
    return k1, k2, k3

def alternating_water_filling(precision,p_w_all,p_p_all,p_w0,p_p0,div,mul,fi_0,cnr_w_d,cnr_p_d,k_w_m,k_p_m):
    p_w = p_w0
    p_p = p_p0
    fi = fi_0
    p_w_new = p_w
    p_p_new = p_p
    while True:
        #print('sum_R',np.sum(np.log2(1+fi)))
        k1,k2,k3=get_approx_parameter(fi, p_w[div], p_p[div], cnr_w_d, cnr_p_d)
        p_w_new, p_p_new=duality_water_filling(div,mul,k1,k2,k3,k_w_m,k_p_m,p_w_all,p_p_all)
        p_w_new=np.clip(p_w_new, a_min=0.001, a_max=None)
        p_p_new=np.clip(p_p_new, a_min=0.001, a_max=None)
        print('p_w_new',np.round(p_w_new,2))
        print('p_w_sum',np.sum(p_w_new))
        print('p_p_new',np.round(p_p_new,2))
        print('p_p_sum',np.sum(p_p_new))
        fi_new = np.exp(k1*np.log(p_w_new[div])+k2*np.log(p_p_new[div])+k3)-1
        fi_new=np.clip(fi_new, a_min=0.00001, a_max=None)
        if np.sqrt(np.linalg.norm(p_w_new-p_w)**2+np.linalg.norm(p_p_new-p_p)**2+np.linalg.norm(fi_new-fi)**2)<precision:
            return p_w_new,p_p_new,fi_new
        else:
            p_w = p_w_new
            p_p = p_p_new
            # print('p_w:',p_w)
            # print('p_p:',p_p)
            fi = fi_new

np.random.seed(11)
facter =1
N =100
sub_band = 20e6 / N
P_w_all = 100
P_p_all = 100
snr_w_ref = np.round(np.clip(np.random.normal(25, 5, size=N), 10, 35), 2)
snr_p_ref = np.round(np.clip(np.random.normal(25, 5, size=N), 10, 35), 2)

cnr_w = np.power(10, snr_w_ref / 10) * N / P_w_all
cnr_p = np.power(10, snr_p_ref / 10) * N / P_p_all

tau_w = gamma.ppf(p_error, m, scale=1 / m)
tau_p = np.exp(s * norm.ppf(p_error) - s ** 2 / 2)

k_w = cnr_w * tau_w
k_p = cnr_p * tau_p

P_w_water = np.clip(water_filling(k_w, P_w_all), a_min=0.0001, a_max=None)
P_p_water = np.clip(water_filling(k_p, P_p_all), a_min=0.0001, a_max=None)

P_w_ave = np.ones(N)*P_w_all/N
P_p_ave = np.ones(N)*P_p_all/N

P_w_o = facter*P_w_water+(1-facter)*P_w_ave
P_p_o = facter*P_p_water+(1-facter)*P_p_ave

plt.figure(0)
plt.bar(np.linspace(1, N, N), P_w_o)
plt.xlabel('subcarriers')
plt.ylabel('p_w')

R_w_mul = np.log2(1+k_w*P_w_water)
R_p_mul = np.log2(1+k_p*P_p_water)

R_w_o = np.log2(1+k_w*P_w_o)
R_p_o = np.log2(1+k_p*P_p_o)

p_sc = p_error_div(R_w_o+R_p_o, P_w_o*cnr_w, P_p_o*cnr_p)
div = np.where(p_sc < p_error)[0]
mul = np.where(p_sc >= p_error)[0]
print('div:', div)
print('mul:', mul)
Nd = len(div)
Nm = len(mul)
k_w_m = k_w[mul]
k_p_m = k_p[mul]
cnr_w_d = cnr_w[div]
cnr_p_d = cnr_p[div]
print('???',1/np.array([]))

p_w0 = np.ones(N)*P_w_all/N
p_p0 = np.ones(N)*P_p_all/N
R_div0 = np.array([outage_R_div(cnr_w_d[i]*p_w0[i], cnr_p_d[i]
                  * p_p0[i], p_error) for i in range(Nd)])
fi_0 = (2**R_div0-1)
precision =0.1
p_w,p_p,fi = alternating_water_filling(precision,P_w_all,P_p_all,p_w0,p_p0,div,mul,fi_0,cnr_w_d,cnr_p_d,k_w_m,k_p_m)
R_old=np.sum(R_w_mul+R_p_mul)*sub_band*1e-6
R_new=(np.sum(np.log2(1+fi))+np.sum(np.log2(1+k_w_m*p_w[mul])+np.log2(1+k_p_m*p_p[mul])))*sub_band*1e-6
print('old',R_old)
print('new',R_new)
plt.figure(1)
plt.bar(np.linspace(1, N, N), snr_w_ref)
plt.xlabel('subcarriers')
plt.ylabel('snr_w')
plt.figure(2)
plt.bar(np.linspace(1, N, N), p_w)
plt.xlabel('subcarriers')
plt.ylabel('p_w')
#plt.show()
