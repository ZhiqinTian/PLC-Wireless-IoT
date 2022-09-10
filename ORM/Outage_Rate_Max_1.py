from scipy.stats import lognorm, gamma, norm
from scipy.optimize import curve_fit
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import numpy as np

p_error = 0.001
s = 1.5
m = 1.3


def p_cdf_hw_2(x):
    return gamma.cdf(x, m, scale=1 / m)


def p_cdf_hp_2(x):
    return lognorm.cdf(x, s, scale=1 / np.exp(s ** 2 / 2))


def func_w(x, a, b):
    return 1-np.exp(-a*x**b)


def func_p(x, a, b):
    return np.exp(-a/x**b)


def p_error_div(R, snr_w, snr_p):
    p_error = p_cdf_hw_2((2**R-1)/snr_w)*p_cdf_hp_2((2**R-1)/snr_p)
    return p_error

def p_error_approx(R,snr_w,a_w,b_w, snr_p,a_p,b_p):
    p_error = func_w((2**R-1)/snr_w,a_w,b_w)*func_p((2**R-1)/snr_p,a_p,b_p)
    return p_error

def p_error_approx1(fi,snr_w,a_w,b_w, snr_p,a_p,b_p):
    p_error = func_w(fi/snr_w,a_w,b_w)*func_p(fi/snr_p,a_p,b_p)
    return p_error

def outage_R_div(snr_w, snr_p, p_error):
    def otage_R_div_eq(log_R):
        return np.log10(p_error_div(2**log_R, snr_w, snr_p))-np.log10(p_error)
    outage_R = np.power(2, fsolve(otage_R_div_eq, [
        np.log2(np.log2(1+max(snr_w, snr_p)))]))
    return outage_R[0]

def outage_R_approx(snr_w,a_w,b_w,snr_p,a_p,b_p, p_error):
    def otage_R_div_eq(log_R):
        return np.log(p_error_approx(2**log_R, snr_w,a_w,b_w, snr_p,a_p,b_p))-np.log(p_error)
    outage_R = np.power(2, fsolve(otage_R_div_eq, [
        np.log2(np.log2(1+max(snr_w, snr_p)))]))
    return outage_R[0]

def outage_R_approx1(snr_w,a_w,b_w,snr_p,a_p,b_p, p_error):
    def otage_R_div_eq(fi_rec):
        return np.log(p_error_approx1(1/fi_rec, snr_w,a_w,b_w, snr_p,a_p,b_p))-np.log(p_error)
    fi = 1/fsolve(otage_R_div_eq, [
        1/max(snr_w, snr_p)])
    return fi[0]

def water_filling(omiga, p_all):
    N = len(omiga)
    flag = np.ones(N)
    p = np.ones(N)
    while True:
        flag[np.where(p < 0)] = 0
        v = np.sum(flag)/(p_all+np.sum((1/omiga)*flag))
        p = (1/v-1/omiga)*flag
        if (p >= 0).all():
            break
    return p


def water_filling_div_mul(omiga_d, omiga_m, p_all):
    Nm = len(omiga_m)
    flag_m = np.ones(Nm)
    p_m = np.ones(Nm)
    while True:
        flag_m[np.where(p_m < 0)] = 0
        v = (np.sum(flag_m)+np.sum(omiga_d))/(p_all+np.sum((1/omiga_m)*flag_m))
        p_d = omiga_d/v
        p_m = (1/v-1/omiga_m)*flag_m
        if (p_m >= 0).all():
            break
    return p_d, p_m


def get_approx_parameter(a_w, b_w, a_p, b_p, fi, p_w, p_p, cnr_w, cnr_p):
    t_w = fi/(p_w*cnr_w)
    t_p = fi/(p_p*cnr_p)
    print('fi',fi)
    print('t_w',t_w)
    print('t_p',t_p)
    e_t_w = np.exp(-a_w*t_w**b_w)
    e_t_p = np.exp(a_p*t_p**(-b_p)+np.log(p_error))
    alpha1 = e_t_w/(e_t_w+e_t_p)
    alpha2 = e_t_p/(e_t_w+e_t_p)
    alpha3 = np.log(e_t_w+e_t_p)-alpha1*(-a_w*t_w**b_w) - \
        alpha2*(a_p*t_p**(-b_p)+np.log(p_error))
    #test=alpha2*(a_p*t_p**(-b_p)+np.log(p_error))-alpha1*(-a_w*t_w**b_w)+alpha3
    beta = t_p**(-b_p)+np.log(p_error)/a_p+alpha3/(alpha2*a_p)
    print('beta',beta)
    beta1 = t_p**(-b_p)/beta
    beta2 = np.log(beta)+beta1*b_p*np.log(t_p)
    omiga = (b_w+b_p*beta1)*(1+fi)
    omiga1 = fi*b_w/omiga
    omiga2 = fi*b_p*beta1/omiga
    omiga3 = omiga1*np.log(cnr_w)+omiga2*np.log(cnr_p)+fi*(beta2+np.log(alpha2*a_p/alpha1/a_w))/omiga+np.log(1+fi)-fi*np.log(fi)/(1+fi)
    return omiga1, omiga2, omiga3

def alternating_water_filling(precision,p_w_all,p_p_all,p_w0,p_p0,div,mul,fi_0,cnr_w_d,cnr_p_d,omiga_w_m,omiga_p_m):
    p_w = p_w0
    p_p = p_p0
    fi = fi_0
    p_w_new = p_w
    p_p_new = p_p
    while True:
        omiga1,omiga2,omiga3=get_approx_parameter(a_w, b_w, a_p, b_p, fi, p_w[div], p_p[div], cnr_w_d, cnr_p_d)
        p_w_new[div],p_w_new[mul] = water_filling_div_mul(omiga1,omiga_w_m,p_w_all)
        p_p_new[div],p_p_new[mul] = water_filling_div_mul(omiga2,omiga_p_m,p_p_all)
        print('p_w_new',p_w_new)
        print('p_w_sum',np.sum(p_w_new))
        print('p_p_new',p_p_new)
        print('p_p_sum',np.sum(p_p_new))
        #p_w_new = np.clip(p_w_new, a_min=0.0001, a_max=None)
        #p_p_new = np.clip(p_p_new, a_min=0.0001, a_max=None)
        print('omiga:',omiga1,omiga2,omiga3)
        fi_new = np.exp(omiga1*np.log(p_w_new[div])+omiga2*np.log(p_p_new[div])+omiga3)-1
        print('p_error_new',func_w(fi_new/(cnr_w_d*p_w_new[div]),a_w,b_w)*func_p(fi_new/(cnr_p_d*p_p_new[div]),a_p,b_p))
        if np.sqrt(np.linalg.norm(p_w_new-p_w)**2+np.linalg.norm(p_p_new-p_p)**2+np.linalg.norm(fi_new-fi)**2)<precision:
            return p_w_new,p_p_new,fi_new
        else:
            p_w = p_w_new
            p_p = p_p_new
            # print('p_w:',p_w)
            # print('p_p:',p_p)
            fi = fi_new
            

plt.figure(1)
xdata1 = [0.6, 0.7, 0.8, 1, 1.3, 1.6, 1.8, 2.4, 2.7]
ydata1 = p_cdf_hw_2(xdata1)
popt_w, pcov = curve_fit(func_w, xdata1, ydata1)
a_w = popt_w[0]
b_w = popt_w[1]
print('w:', a_w, b_w)
xdata2 = [0.01, 0.02, 0.05, 0.08, 0.1, 0.3, 0.4, 0.5,
          0.6, 0.7, 0.8, 1, 1.3, 1.6, 1.8, 2.4, 2.7, 3.5]
ydata2 = p_cdf_hp_2(xdata2)
popt_p, pcov = curve_fit(func_p, xdata2, ydata2)
a_p = popt_p[0]
b_p = popt_p[1]

# print('p:', a_p, b_p)
# x = np.linspace(0.05, 10, 100)
# # plt.plot(x,np.exp(-5*x)+np.exp(-0.2/x))
# plt.plot(x, np.log(p_cdf_hw_2(x)), label='p_w')
# plt.plot(x, np.log(1-np.exp(-a_w*x**b_w)), 'b--')
# plt.legend()
# plt.figure(2)
# plt.plot(x, -a_p/x**b_p, 'g--')
# plt.plot(x, np.log(p_cdf_hp_2(x)), label='p_p')
# plt.legend()
# plt.show()
np.random.seed(1)
N =20
sub_band = 20e6 / N
P_w_all = 100
P_p_all = 100
snr_w_ref = np.round(np.clip(np.random.normal(20, 10, size=N), 10, 40), 2)
snr_p_ref = np.round(np.clip(np.random.normal(20, 10, size=N), 10, 40), 2)
plt.plot(np.linspace(1, N, N), snr_p_ref)

cnr_w = np.power(10, snr_w_ref / 10) * N / P_w_all
cnr_p = np.power(10, snr_p_ref / 10) * N / P_p_all

k_p = np.exp(s * norm.ppf(p_error) - s ** 2 / 2)
k_w = gamma.ppf(p_error, m, scale=1 / m)
omiga_w = cnr_w * k_w
omiga_p = cnr_p * k_p
P_w = np.clip(water_filling(omiga_w, P_w_all), a_min=0.0001, a_max=None)
P_p = np.clip(water_filling(omiga_p, P_p_all), a_min=0.0001, a_max=None)
R_w = np.log2(1+omiga_w*P_w)
R_p = np.log2(1+omiga_p*P_p)
R_w_all = np.sum(R_w)
R_p_all = np.sum(R_p)
p_sc = p_error_div(R_w+R_p, P_w*cnr_w, P_p*cnr_p)
div = np.where(p_sc < p_error)[0]
mul = np.where(p_sc >= p_error)[0]
print('div:', div)
print('mul:', mul)
Nd = len(div)
Nm = len(mul)
omiga_w_m = omiga_w[mul]
omiga_p_m = omiga_p[mul]
cnr_w_d = cnr_w[div]
cnr_p_d = cnr_p[div]

p_w0 = np.ones(N)*P_w_all/N
p_p0 = np.ones(N)*P_p_all/N
R_div0 = np.array([outage_R_div(cnr_w_d[i]*p_w0[i], cnr_p_d[i]
                  * p_p0[i], p_error) for i in range(Nd)])
# R_div0_ap = np.array([outage_R_approx(cnr_w_d[i]*p_w0[i],a_w,b_w, cnr_p_d[i]
#                   * p_p0[i],a_p,b_p, p_error) for i in range(Nd)])
fi_ap = 0.9*np.array([outage_R_approx1(cnr_w_d[i]*p_w0[i],a_w,b_w, cnr_p_d[i]
                  * p_p0[i],a_p,b_p, p_error) for i in range(Nd)])
R_div0_ap =np.log2(1+fi_ap)
p_error = p_error_div(R_div0_ap,cnr_w_d*p_w0[div],cnr_p_d*p_p0[div])
p_error_ap = func_w(fi_ap/(cnr_w_d*p_w0[div]),a_w,b_w)*func_p(fi_ap/(cnr_p_d*p_p0[div]),a_p,b_p)
print('p_error',p_error)
print('p_error_ap',p_error_ap)
print('R',R_div0)
print('R_ap:',R_div0_ap)
fi_0 = fi_ap
precision =0.1
# fi=np.linspace(10,100,1000)
# plt.figure(7)
# test=np.log(p_error_approx1(fi, 100,a_w,b_w, 100,a_p,b_p))
# plt.plot(1/fi,test)
# plt.show()
#print(fi_0)
#omiga1,omiga2,omiga3=get_approx_parameter(a_w, b_w, a_p, b_p, fi_0, p_w0[div], p_p0[div], cnr_w_d, cnr_p_d)
#fi = np.exp(omiga1*np.log(p_w0[div])+omiga2*np.log(p_p0[div])+omiga3)-1
p_w,p_p,fi = alternating_water_filling(precision,P_w_all,P_p_all,p_w0,p_p0,div,mul,fi_0,cnr_w_d,cnr_p_d,omiga_w_m,omiga_p_m)
print('p_w',np.round(p_w,2))
print('p_p',np.round(p_p,2))
# omiga_d=np.array([1,2,3])
# omiga_m=np.array([0.01,3,10,4,6])
# Nd=len(omiga_d)
# Nm=len(omiga_m)
# p_all=100
# p_d,p_m=water_filling(omiga_d,omiga_m,p_all)
# print(p_d,p_m)
# plt.figure(1)
# plt.bar(np.linspace(1, Nd,Nd ), p_d)
# plt.xlabel('subcarriers')
# plt.ylabel('pd')
plt.figure(1)
plt.bar(np.linspace(1, N, N), snr_w_ref)
plt.xlabel('subcarriers')
plt.ylabel('snr_w')
plt.show()
