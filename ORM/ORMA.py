from tabnanny import verbose
from scipy.stats import lognorm, gamma, norm,ncx2
from scipy.optimize import fsolve
import numpy as np
import cvxpy as cp
import numpy as np
import matplotlib.pyplot as plt

class ORMA:
    def __init__(self, p_error=0.001,e_2=0.05, sm=None):
        self.p_error = p_error
        self.sm = sm
        self.e_2=e_2
    
    def set_pe(self,pe):
        self.p_error=pe

    def Outage_Rate_Max_Algorithm(self, cnr_w, hf_2, cnr_p, s_p, P_w_all, P_p_all):
        div, mul = self.AADS(cnr_w, hf_2, cnr_p, s_p, P_w_all, P_p_all)
        p_w, p_p, R_all , iter_num= self.optimal_power_allocation(
            cnr_w, hf_2, cnr_p, s_p, P_w_all, P_p_all, div, mul)
        # print('div',div)
        # print('mul',mul)
        # print('R_all',R_all)
        return {'div':div, 'mul':mul,'p_w': p_w,'p_p': p_p,'R_all': R_all,'iter_num':iter_num}
    
    # def Outage_Rate_Only_Mul(self, cnr_w, hf_2, cnr_p, s_p, P_w_all, P_p_all):
    #     p_w=self.water_filling(cnr_w*hf_2, P_w_all)
    #     p_p=self.water_filling(cnr_p,P_p_all)
    #     R_w=np.sum(np.log2(1+self.tau_w*cnr_w*p_w))*self.sm.spacing*1e-6
    #     R_p=np.sum(np.log2(1+self.tau_p*cnr_p*p_p))*self.sm.spacing*1e-6
    #     R_all=R_w+R_p
    #     return {'p_w': p_w,'p_p': p_p,'R_all': R_all}
    
    # def Outage_Rate_Only_Div(self, cnr_w, hf_2, cnr_p, s_p, P_w_all, P_p_all):
    #     div, mul = list(range(len(cnr_w))),[]
    #     p_w, p_p, R_all , iter_num= self.optimal_power_allocation(
    #         cnr_w, hf_2, cnr_p, s_p, P_w_all, P_p_all, div, mul)
    #     return {'p_w': p_w,'p_p': p_p,'R_all': R_all,'iter_num':iter_num}

    def pdf_hw_2(self,x,hf_2):
        return ncx2.pdf(x, 2, hf_2*2/self.e_2, scale=self.e_2/2)

    def cdf_hw_2(self,x,hf_2):
        return ncx2.cdf(x, 2, hf_2*2/self.e_2, scale=self.e_2/2)

    def ppf_hw_2(self,x,hf_2):
        return ncx2.ppf(x, 2, hf_2*2/self.e_2, scale=self.e_2/2)

    def pdf_hp_2(self,x,s):
        return lognorm.pdf(x, s, scale=1 / np.exp(s ** 2 / 2))

    def cdf_hp_2(self,x,s):
        return lognorm.cdf(x, s, scale=1 / np.exp(s ** 2 / 2))

    def ppf_hp_2(self,x,s):
        return lognorm.ppf(x, s, scale=1 / np.exp(s ** 2 / 2))

    def p_error_div(self,R, snr_w, hf_2, snr_p,s):
        p_error = self.cdf_hw_2((2**R-1)/snr_w,hf_2)*self.cdf_hp_2((2**R-1)/snr_p,s)
        return p_error

    def outage_R_div(self, p_outage, snr_w, hf_2, snr_p, s):
      
        C_w=np.log2(1+hf_2*snr_w)
        C_p=np.log2(1+snr_p)
        # print('rw',snr_w)
        # print('rp',snr_p)
        Ra=0
        Rb=C_w+C_p
        dealta_log_p=1
        while np.abs(dealta_log_p)>0.05:
            Rc=(Ra+Rb)/2
            dealta_log_p=np.log(self.p_error_div(Rc, snr_w, hf_2, snr_p,s))-np.log(p_outage)
            if dealta_log_p<0:
                Ra=Rc
            else:
                Rb=Rc
        return Rc
    
    def outage_R_mul(self, p_outage, snr_w, hf_2, snr_p, s):
        k1=self.ppf_hw_2(p_outage, hf_2)*self.pdf_hw_2(self.ppf_hw_2(p_outage, hf_2), hf_2)/p_outage
        k2=self.ppf_hp_2(p_outage, s)*self.pdf_hp_2(self.ppf_hp_2(p_outage, s), s)/p_outage
        
        kw=self.ppf_hw_2(k2/(k1+k2)*p_outage, hf_2)
        kp=self.ppf_hp_2(k1/(k1+k2)*p_outage, s)
        
        R_w = np.log2(1+kw*snr_w)
        R_p = np.log2(1+kp*snr_p)
        return R_w+R_p

    def water_filling(self, omiga, p_all):
        N = len(omiga)
        flag = np.ones(N)
        p = np.ones(N)
        omiga = np.clip(omiga, a_min=0.00001, a_max=None)
        while True:
            flag[np.where(p < 0)] = 0
            miu = np.sum(flag)/(p_all+np.sum((1/omiga)*flag))
            p = (1/miu-1/omiga)*flag
            if (p >= 0).all():
                break
        return np.abs(p)

    def water_filling_div_mul(self, omiga_d, omiga_m, p_all):
        Nm = len(omiga_m)
        flag_m = np.ones(Nm)
        p_m = np.ones(Nm)
        omiga_d = np.clip(omiga_d, a_min=0.00001, a_max=None)
        omiga_m = np.clip(omiga_m, a_min=0.00001, a_max=None)
        t=1
        while True:
            t=t+1
            flag_m[np.where(p_m < 0)] = 0
            miu = (np.sum(flag_m)+np.sum(omiga_d)) / \
                (p_all+np.sum((1/omiga_m)*flag_m))
            p_d = omiga_d/miu
            p_m = (1/miu-1/omiga_m)*flag_m
            
            if (p_m >= 0).all():
                break
        return p_d, np.abs(p_m)
    
    def AADS(self, cnr_w, hf_2, cnr_p, s_p, P_w_all, P_p_all):
        N = self.sm.num_subcarriers
        snr_w=np.ones(N)*P_w_all/N*cnr_w
        snr_p=np.ones(N)*P_p_all/N*cnr_p
        R_mul_o = self.outage_R_mul(self.p_error,snr_w,hf_2,snr_p,s_p)
        #print('R_mul_o',R_mul_o)
        p_sc = self.p_error_div(R_mul_o, snr_w, hf_2, snr_p,s_p)
        div = np.where(p_sc < self.p_error)[0]
        mul = np.where(p_sc >= self.p_error)[0]
        return div, mul

    def subopt_approx_parameter(self, fi, p_w_d, p_p_d, cnr_w, hf_2_d, cnr_p_d, s_d):
        x_w = fi/(p_w_d*cnr_w)
        x_p = fi/(p_p_d*cnr_p_d)
        fw = self.pdf_hw_2(x_w,hf_2_d)
        Fw = np.clip(self.cdf_hw_2(x_w, hf_2_d), a_min=1e-8, a_max=None)
        fp = self.pdf_hp_2(x_p, s_d)
        Fp = np.clip(self.cdf_hp_2(x_p, s_d), a_min=1e-8, a_max=None)
        alpha1 = fw*x_w/Fw
        alpha2 = fp*x_p/Fp
        beta1 = np.log(Fw)-alpha1*np.log(x_w)
        beta2 = np.log(Fp)-alpha2*np.log(x_p)
        omiga1 = alpha1/(alpha1+alpha2)
        omiga2 = alpha2/(alpha1+alpha2)
        omiga3 = omiga1*np.log(cnr_w)+omiga2*np.log(cnr_p_d)-(beta1+beta2) / \
            (alpha1+alpha2)+np.log(self.p_error)/(alpha1+alpha2)
        k0 = fi/(fi+1)
        k1 = np.clip(omiga1*k0, a_min=1e-8, a_max=None)
        k2 = np.clip(omiga2*k0, a_min=1e-8, a_max=None)
        k3 = np.clip(omiga3*k0+np.log(1+fi)-k0*np.log(fi), a_min=1e-8, a_max=None)
        return k1, k2, k3
    
    def get_approx_parameter(self, fi, p_w_d, p_p_d, cnr_w, hf_2_d, cnr_p_d, s_d):
        x_w = fi/(p_w_d*cnr_w)
        x_p = fi/(p_p_d*cnr_p_d)
        fw = self.pdf_hw_2(x_w,hf_2_d)
        Fw=self.cdf_hw_2(x_w, hf_2_d)
        # print('Fw',Fw)
        Fw = np.clip(Fw, a_min=1e-8, a_max=None)
        fp = self.pdf_hp_2(x_p, s_d)
        Fp=self.cdf_hp_2(x_p, s_d)
        # print('Fp',Fp)
        Fp = np.clip(Fp, a_min=1e-8, a_max=None)
        alpha1 = fw/Fw
        alpha2 = fp/Fp
        beta1 = np.log(Fw)-alpha1*x_w
        beta2 = np.log(Fp)-alpha2*x_p
        # print('beta1',beta1)
        # print('beta2',beta2)
        omega1 = alpha1/cnr_w
        omega2 = alpha2/cnr_p_d
        omega3 = (beta1+beta2)-np.log(self.p_error)
        # print('w3',omega3)
        omega1 = np.clip(omega1, a_min=1e-30, a_max=None)
        omega2 = np.clip(omega2, a_min=1e-30, a_max=None)
        #omega3 = np.clip(omega3, a_min=1e-10, a_max=None)
        
        k_d = fi/(fi+1)
        Omega_d = np.log(1+fi)-k_d*np.log(fi)
        #Omega_d = np.clip(Omega_d, a_min=1e-8, a_max=None)
        # print('k_d',k_d)
        # print('Omega_d',Omega_d)
        # print('omega1',omega1)
        # print('omega2',omega2)
        # print('omega3',omega3)
        return k_d,Omega_d,omega1,omega2,omega3
    
    
    # def SICAA(self, precision, p_w_all, p_p_all, p_w0, p_p0, div, mul, fi_0, cnr_w, hf_2_d, cnr_p_d, s_d, k_w_m, k_p_m):
    #     p_w = p_w0
    #     p_p = p_p0
    #     fi = fi_0
    #     p_w_new = p_w
    #     p_p_new = p_p
    #     t=0
    #     IWFA_t=[]
    #     while t<100:
    #         t+=1
    #         #print('tt:',t)
    #         k1, k2, k3 = self.get_approx_parameter(
    #             fi, p_w[div], p_p[div], cnr_w, hf_2_d, cnr_p_d, s_d)
    #         p_w_new, p_p_new,iwfa_t = self.IWFA(
    #             div, mul, k1, k2, k3, k_w_m, k_p_m, p_w_all, p_p_all)
    #         #print('twf:',t)
    #         IWFA_t.append(iwfa_t)
    #         p_w_new = np.clip(p_w_new, a_min=0.001, a_max=None)
    #         p_p_new = np.clip(p_p_new, a_min=0.001, a_max=None)
    #         fi_new = np.exp(
    #             k1*np.log(p_w_new[div])+k2*np.log(p_p_new[div])+k3)-1
    #         fi_new = np.clip(fi_new, a_min=0.00001, a_max=None)
    #         if np.sqrt(np.linalg.norm(p_w_new-p_w)**2+np.linalg.norm(p_p_new-p_p)**2+np.linalg.norm(fi_new-fi)**2) < precision:
    #             print('SICAA_t:',t)
    #             return p_w_new, p_p_new, fi_new,{'SICAA_t':t,'IWFA_t':IWFA_t}
    #         else:
    #             p_w = p_w_new
    #             p_p = p_p_new
    #             fi = fi_new
    #     print('SICAA_t:',t)
    #     return p_w_new, p_p_new, fi_new,{'SICAA_t':t,'IWFA_t':IWFA_t}
    
    def SICAA(self, precision, p_w_all, p_p_all, p_w0, p_p0, div, mul, fi_0, cnr_w, hf_2_d, cnr_p_d, s_d, k_w_m, k_p_m):
        p_w = p_w0
        p_p = p_p0
        fi = fi_0
        sub_band=self.sm.spacing
        R=(np.sum(np.log2(1+fi_0))+np.sum(np.log2(1+k_w_m *
                                                      p_w0[mul])+np.log2(1+k_p_m*p_p0[mul])))*sub_band*1e-6
        p_w_new = p_w
        p_p_new = p_p
        Nd=len(div)
        Nm=len(mul)
        t=0
        # print('k_w_m',k_w_m)
        # print('k_p_m', k_p_m)
        while t<10:
            t+=1
            #print('tt:',t)
            k_d,_,omega1,omega2,omega3 = self.get_approx_parameter(
                fi, p_w[div], p_p[div], cnr_w, hf_2_d, cnr_p_d, s_d)
            y_w_d,y_p_d,z_d,p_w_m,p_p_m=self.convex_slution(k_d,k_w_m, k_p_m,omega1,omega2,omega3, p_w_all, p_p_all,Nd,Nm)
            # print(y_w_d)
            p_w_new[div] = np.exp(y_w_d)
            p_p_new[div] = np.exp(y_p_d)
            p_w_new[mul] = p_w_m
            p_p_new[mul] = p_p_m
            fi_new = np.exp(z_d)
            p_w_new = np.clip(p_w_new, a_min=0.001, a_max=None)
            p_p_new = np.clip(p_p_new, a_min=0.001, a_max=None)
            fi_new = np.clip(fi_new, a_min=0.00001, a_max=None)
            R_new = (np.sum(np.log2(1+fi_new))+np.sum(np.log2(1+k_w_m *
                                                      p_w_new[mul])+np.log2(1+k_p_m*p_p_new[mul])))*sub_band*1e-6
            #print('R',R_new)
            # if np.sqrt(np.linalg.norm(p_w_new-p_w)**2+np.linalg.norm(p_p_new-p_p)**2+np.linalg.norm(fi_new-fi)**2) < precision:
            if np.linalg.norm(R_new-R) < precision:
                print('SICAA_t:',t)
                return p_w_new, p_p_new, fi_new,{'SICAA_t':t}
            else:
                p_w = p_w_new
                p_p = p_p_new
                fi = fi_new
                R=R_new
        print('SICAA_t:',t)
        return p_w_new, p_p_new, fi_new,{'SICAA_t':t}
    
    def convex_slution(self,k_d,k_w_m,k_p_m,w1_d,w2_d,w3_d,Pw,Pp,Nd,Nm):
        
        mosek_para = {"MSK_DPAR_INTPNT_CO_TOL_REL_GAP": 1e-1}
        #print(Nd,Nm)
        if Nd and Nm:
            # print('k_d',k_d)
            # print('k_w_m',k_w_m)
            # print('k_p_m',k_p_m)
            # print('w1_d',w1_d)
            # print('w2_d',w2_d)
            # print('w3_d',w3_d)
            y_w_d = cp.Variable(Nd)
            p_w_m = cp.Variable(Nm)
            y_p_d = cp.Variable(Nd)
            p_p_m = cp.Variable(Nm)
            z_d= cp.Variable(Nd)
            objective = cp.Maximize(cp.sum(cp.multiply(k_d,z_d))+cp.sum(cp.log(1+cp.multiply(k_w_m,p_w_m))+cp.log(1+cp.multiply(k_p_m,p_p_m))))
            # constraints = [ cp.multiply(w1_d,cp.exp(z_d-y_w_d))+cp.multiply(w2_d,cp.exp(z_d-y_p_d))+w3_d<=0,
            #             cp.sum(cp.exp(y_w_d))+cp.sum(p_w_m)<= Pw,
            #             cp.sum(cp.exp(y_p_d))+cp.sum(p_p_m)<= Pp,
            #             -cp.multiply(k_d,z_d)-omega_d<=0,
            #             -p_w_m<=0,
            #             -p_p_m<=0]
            constraints = [ cp.multiply(w1_d,cp.exp(z_d-y_w_d))+cp.multiply(w2_d,cp.exp(z_d-y_p_d))+w3_d<=0,
                        cp.sum(cp.exp(y_w_d))+cp.sum(p_w_m)<= Pw,
                        cp.sum(cp.exp(y_p_d))+cp.sum(p_p_m)<= Pp,
                        -p_w_m<=0,
                        -p_p_m<=0]
            prob = cp.Problem(objective, constraints)
            prob.solve(solver=cp.MOSEK)
            #prob.solve()
            #print('yd_value',y_w_d.value)
            return y_w_d.value,y_p_d.value,z_d.value,p_w_m.value,p_p_m.value
        elif not Nd:
            p_w_m_new=self.water_filling(k_w_m,Pw)
            p_p_m_new=self.water_filling(k_p_m,Pp)
            return [],[],[],p_w_m_new,p_p_m_new
        else:
            y_w_d = cp.Variable(Nd)
            y_p_d = cp.Variable(Nd)
            z_d= cp.Variable(Nd)
            # print('dim',w2_d.shape)
            # t=cp.multiply(w2_d,cp.exp(z_d-y_p_d))
            objective = cp.Maximize(cp.sum(cp.multiply(k_d,z_d)))
            # constraints = [ cp.multiply(w1_d,cp.exp(z_d-y_w_d))+cp.multiply(w2_d,cp.exp(z_d-y_p_d))+w3_d<=0,
            #             cp.sum(cp.exp(y_w_d))<= Pw,
            #             cp.sum(cp.exp(y_p_d))<= Pp,
            #             -cp.multiply(k_d,z_d)-omega_d<=0]
            constraints = [ cp.multiply(w1_d,cp.exp(z_d-y_w_d))+cp.multiply(w2_d,cp.exp(z_d-y_p_d))+w3_d<=0,
                        cp.sum(cp.exp(y_w_d))<= Pw,
                        cp.sum(cp.exp(y_p_d))<= Pp,]
            prob = cp.Problem(objective, constraints)
            #prob.solve(solver=cp.MOSEK,mosek_params=mosek_para)
            prob.solve(solver=cp.MOSEK)
            #print('Nd',Nd)
            return y_w_d.value,y_p_d.value,z_d.value,[],[]
    
    def IWFA(self, div, mul, k1, k2, k3, k_w_m, k_p_m, p_w_all, p_p_all):
        Nd = len(div)
        Nm = len(mul)
        p_w_new = np.ones(Nd+Nm)
        p_p_new = np.ones(Nd+Nm)
        v_d = 0.01*np.ones(Nd)
        k_v = 0.01
        t=0
        while t<50:
            t+=1
            print('wftt',t)
            p_w_new[div], p_w_new[mul] = self.water_filling_div_mul(
                k1+v_d, k_w_m, p_w_all)
            p_p_new[div], p_p_new[mul] = self.water_filling_div_mul(
                k2+v_d, k_p_m, p_p_all)
            p_w_new = np.clip(p_w_new, a_min=0.001, a_max=None)
            p_p_new = np.clip(p_p_new, a_min=0.001, a_max=None)
            v_d_new = np.clip(
                v_d-k_v*(k1*np.log(p_w_new[div])+k2*np.log(p_p_new[div])+k3), a_min=0, a_max=None)
            if np.linalg.norm(v_d_new-v_d) < 0.005:
                #print('IWFA_t:',t)
                return p_w_new, p_p_new,t
            else:
                v_d = v_d_new
        #print('IWFA_t:',t)
        return p_w_new, p_p_new,t
    
    def optimal_power_allocation(self, cnr_w,hf_2, cnr_p, s_p, P_w_all, P_p_all, div, mul):
        N = self.sm.num_subcarriers
        sub_band=self.sm.spacing
        Nd = len(div)
        k1=self.ppf_hw_2(self.p_error, hf_2)*self.pdf_hw_2(self.ppf_hw_2(self.p_error, hf_2), hf_2)/self.p_error
        k2=self.ppf_hp_2(self.p_error, s_p)*self.pdf_hp_2(self.ppf_hp_2(self.p_error, s_p), s_p)/self.p_error
        kw=self.ppf_hw_2(k2/(k1+k2)*self.p_error, hf_2)
        kp=self.ppf_hp_2(k1/(k1+k2)*self.p_error, s_p)
        k_w=kw*cnr_w
        k_p=kp*cnr_p
        k_w_m = k_w[mul]
        k_p_m = k_p[mul]
        cnr_p_d = cnr_p[div]
        hf_2_d=hf_2[div]
        s_p_d=s_p[div]
        p_w0 = np.ones(N)*P_w_all/N
        p_p0 = np.ones(N)*P_p_all/N
        R_div0 = np.array([self.outage_R_div(self.p_error, cnr_w*p_w0[i], hf_2_d[i], cnr_p_d[i] * p_p0[i],s_p_d[i]) for i in range(Nd)])
        # print('cnr_w',cnr_w)
        # print('cnr_p_d',cnr_p_d)
        # print('hf_2_d',hf_2_d)
        # print('s_p',hf_2_d)
        fi_0 = (2**R_div0-1)
        # print('fi_0',fi_0)
        precision = 0.1
        p_w, p_p, fi,iter_num = self.SICAA(
            precision, P_w_all, P_p_all, p_w0, p_p0, div, mul, fi_0, cnr_w, hf_2_d, cnr_p_d, s_p_d, k_w_m, k_p_m)
        R_new = (np.sum(np.log2(1+fi))+np.sum(np.log2(1+k_w_m *
                                                      p_w[mul])+np.log2(1+k_p_m*p_p[mul])))*sub_band*1e-6
        return p_w, p_p, R_new,iter_num
    
    def suboptimal_power_allocation(self, cnr_w,hf_2, cnr_p, s_p, P_w_all, P_p_all, div, mul):
        N = self.sm.num_subcarriers
        #sub_band=self.sm.spacing
        Nd = len(div)
        k1=self.ppf_hw_2(self.p_error, hf_2)*self.pdf_hw_2(self.ppf_hw_2(self.p_error, hf_2), hf_2)/self.p_error
        k2=self.ppf_hp_2(self.p_error, s_p)*self.pdf_hp_2(self.ppf_hp_2(self.p_error, s_p), s_p)/self.p_error
        kw=self.ppf_hw_2(k2/(k1+k2)*self.p_error, hf_2)
        kp=self.ppf_hp_2(k1/(k1+k2)*self.p_error, s_p)
        k_w=kw*cnr_w
        k_p=kp*cnr_p
        k_w_m = k_w[mul]
        k_p_m = k_p[mul]
        cnr_p_d = cnr_p[div]
        hf_2_d=hf_2[div]
        s_p_d=s_p[div]
        p_w0 = np.ones(N)*P_w_all/N
        p_p0 = np.ones(N)*P_p_all/N
        R_div0 = np.array([self.outage_R_div(self.p_error, cnr_w*p_w0[i], hf_2_d[i], cnr_p_d[i] * p_p0[i],s_p_d[i]) for i in range(Nd)])
        fi_0 = (2**R_div0-1)
        k1, k2, k3 = self.subopt_approx_parameter(fi_0, p_w0[div], p_p0[div], cnr_w, hf_2_d, cnr_p_d, s_p_d)
        p_w=p_w0
        p_p=p_p0
        p_w[div], p_w[mul] = self.water_filling_div_mul(k1, k_w_m, P_w_all)
        p_p[div], p_p[mul] = self.water_filling_div_mul(k2, k_p_m, P_p_all)
        return p_w, p_p
    
    def uniform_power_allocation(self, P_w_all, P_p_all):
        N = self.sm.num_subcarriers
        p_w0 = np.ones(N)*P_w_all/N
        p_p0 = np.ones(N)*P_p_all/N
        p_w, p_p  = p_w0, p_p0
        return p_w, p_p
    
    def water_filling_allocation(self, cnr_w, hf_2, cnr_p, P_w_all, P_p_all):
        p_w = self.water_filling(cnr_w*hf_2,P_w_all)
        p_p = self.water_filling(cnr_p,P_p_all)
        return p_w, p_p
    
    def single_link(self, cnr_w, hf_2, cnr_p,s_p, p_w, p_p):
        sub_band=self.sm.spacing
        kw=self.ppf_hw_2(self.p_error, hf_2)
        kp=self.ppf_hp_2(self.p_error, s_p)
        R_w=np.log2(1+kw*cnr_w*p_w)*sub_band*1e-6
        R_p=np.log2(1+kp*cnr_p*p_p)*sub_band*1e-6
        return R_w, R_p
        
    def outage_rate(self, cnr_w, hf_2, cnr_p, s_p, p_w, p_p, div, mul):
        sub_band=self.sm.spacing
        k1=self.ppf_hw_2(self.p_error, hf_2)*self.pdf_hw_2(self.ppf_hw_2(self.p_error, hf_2), hf_2)/self.p_error
        k2=self.ppf_hp_2(self.p_error, s_p)*self.pdf_hp_2(self.ppf_hp_2(self.p_error, s_p), s_p)/self.p_error
        kw=self.ppf_hw_2(k2/(k1+k2)*self.p_error, hf_2)
        kp=self.ppf_hp_2(k1/(k1+k2)*self.p_error, s_p)
        k_w=kw*cnr_w
        k_p=kp*cnr_p
        k_w_m = k_w[mul]
        k_p_m = k_p[mul]
        #print('ok1')
        R_div = np.sum(np.array([self.outage_R_div(self.p_error, cnr_w*p_w[i], hf_2[i], cnr_p[i]*p_p[i], s_p[i]) for i in div]))*sub_band*1e-6
        #print('ok2')
        R_w_m = np.sum(np.log2(1+k_w_m *p_w[mul]))*sub_band*1e-6
        R_p_m = np.sum(np.log2(1+k_p_m *p_p[mul]))*sub_band*1e-6
        R_outage = (R_div+R_w_m+R_p_m)
        return R_outage