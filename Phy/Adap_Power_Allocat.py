import numpy as np
from scipy.stats import lognorm, gamma, norm
from scipy.optimize import fsolve


class Adap_PA:
    def __init__(self, m=1.2, s=1.5, p_error=1e-3, sm=None):
        self.m = m
        self.s = s
        self.p_error = p_error
        self.sm = sm

    def set(self, s, m, p_error, sm):
        self.s = s
        self.m = m
        self.p_error = p_error
        self.sm = sm

    def p_pdf_hw_2(self, x):
        return gamma.pdf(x, self.m, scale=1 / self.m)

    def p_pdf_hp_2(self, x):
        return lognorm.pdf(x, self.s, scale=1 / np.exp(self.s ** 2 / 2))

    def p_cdf_hw_2(self, x):
        return gamma.cdf(x, self.m, scale=1 / self.m)

    def p_cdf_hp_2(self, x):
        return lognorm.cdf(x, self.s, scale=1 / np.exp(self.s ** 2 / 2))

    def p_outage_sc_dive(self, R, snr_w, m, snr_p, s):
        def p_outage_wireless(R, snr_w, m):
            return gamma.cdf((2 ** R - 1) / snr_w, m, scale=1 / m)

        def p_outage_plc(R, snr_p, s):
            return lognorm.cdf((2 ** R - 1) / snr_p, s, scale=1 / np.exp(s ** 2 / 2))

        p_outage = p_outage_wireless(R, snr_w, m) * p_outage_plc(R, snr_p, s)
        return p_outage

    def R_sc(self, p_error, snr_w, m, snr_p, s):

        def sc_dive_R_max(log_R):
            return np.log10(self.p_outage_sc_dive(2**log_R, snr_w, m, snr_p, s))-np.log10(p_error)
        R_sc = np.power(2, fsolve(sc_dive_R_max, [
            np.log2(np.log2(1+max(snr_w, snr_p)))]))
        return R_sc[0]

    def adap_PA(self, cnr_w, tx_power_w, cnr_p, tx_power_p):
        N = self.sm.num_subcarriers
        sub_band = self.sm.spacing
        P_w_all = tx_power_w
        P_p_all = tx_power_p

        k_p = np.exp(self.s * norm.ppf(self.p_error) - self.s ** 2 / 2)
        k_w = gamma.ppf(self.p_error, self.m, scale=1 / self.m)
        factor_w = cnr_w * k_w
        factor_p = cnr_p * k_p
        P_w = np.clip((P_w_all + np.sum(1 / factor_w)) /
                      N - 1 / factor_w, a_min=1e-6, a_max=None)
        P_p = np.clip((P_p_all + np.sum(1 / factor_p)) / N -
                      1 / factor_p, a_min=1e-6, a_max=None)
        P_w = P_w * P_w_all / np.sum(P_w)
        P_p = P_p * P_p_all / np.sum(P_p)
        R_w = np.log2(1+factor_w*P_w)
        R_p = np.log2(1+factor_p*P_p)
        R_w_all = np.sum(R_w)
        R_p_all = np.sum(R_p)
        R_reuse = (R_w_all+R_p_all)*sub_band

        p_sc = self.p_outage_sc_dive(
            R_w+R_p, P_w*cnr_w, self.m, P_p*cnr_p, self.s)
        d = np.where(p_sc < self.p_error)[0]
        r = np.where(p_sc >= self.p_error)[0]

        Nd = len(d)
        Nr = len(r)
        factor_w_r = factor_w[r]
        factor_p_r = factor_p[r]
        cnr_w_d = cnr_w[d]
        cnr_p_d = cnr_p[d]

        def div_reuse_equation(x):
            X_w_d = x[0:Nd]
            X_p_d = x[Nd:2*Nd]
            P_w_d = np.exp(X_w_d)
            P_p_d = np.exp(X_p_d)
            X_fi = x[2*Nd:3*Nd]

            snr_div = np.exp(X_fi)
            snr_w_d = P_w_d*cnr_w_d
            snr_p_d = P_p_d*cnr_p_d

            fw = self.p_pdf_hw_2(snr_div/snr_w_d)
            fp = self.p_pdf_hp_2(snr_div/snr_p_d)
            Fw = self.p_cdf_hw_2(snr_div/snr_w_d)
            Fp = self.p_cdf_hp_2(snr_div/snr_p_d)

            delta_w = np.sum(P_w_d)-P_w_all-np.sum(1/factor_w_r)
            delta_p = np.sum(P_p_d)-P_p_all-np.sum(1/factor_p_r)
            yw_part = (1+snr_div)*(1+snr_w_d * fp*Fw/(snr_p_d*(fw+0.01)*Fp))
            yw = P_w_d*yw_part/snr_div+delta_w/Nr+0.3*X_w_d+0.5*X_fi
            yp_part = (1+snr_div)*(1+snr_p_d * fw*Fp/(snr_w_d*(fp+0.01)*Fw))
            yp = P_p_d*yp_part/snr_div+delta_p/Nr+0.3*X_p_d+0.5*X_fi

            log_pe = np.log(self.p_error)-np.log(1e-12+Fp)-np.log(1e-12+Fw)

            result = np.concatenate((yw, yp, log_pe), axis=0)
            return result

        if d.size > 0 and r.size > 0:
            x_w = np.log(np.ones(Nd)*P_w_all/N)
            x_p = np.log(np.ones(Nd)*P_p_all/N)
            p_w_0 = np.exp(x_w)
            p_p_0 = np.exp(x_p)
            R_W0 = np.log2(1+p_w_0*factor_w[d])
            R_P0 = np.log2(1+p_p_0*factor_p[d])
            x_fi = np.log(np.power(2, R_W0+R_P0)-1)
            x0 = np.concatenate(
                (x_w, x_p, x_fi), axis=0)
            result = np.round(fsolve(div_reuse_equation, x0), 2)

            P_w_d = np.exp(result[0:Nd])
            P_p_d = np.exp(result[Nd:2*Nd])
            #snr_div = np.exp(result[2*Nd:3*Nd])
            if r.size > 0:
                P_w_r = (P_w_all-np.sum(P_w_d)+np.sum(1/factor_w_r)) / \
                    Nr - 1/factor_w_r
                P_p_r = (P_p_all-np.sum(P_p_d)+np.sum(1/factor_p_r)) / \
                    Nr - 1/factor_p_r

            P_w_new = np.zeros(N)
            P_p_new = np.zeros(N)

            P_w_new[d] = P_w_d
            P_p_new[d] = P_p_d
            if r.size > 0:
                P_w_new[r] = P_w_r
                P_p_new[r] = P_p_r

            P_w_new = np.clip(P_w_new, a_min=1e-6, a_max=None)
            P_p_new = np.clip(P_p_new, a_min=1e-6, a_max=None)
            P_w_new = P_w_new * P_w_all / np.sum(P_w_new)
            P_p_new = P_p_new * P_p_all / np.sum(P_p_new)

            P_w_d_new = P_w_new[d]
            P_p_d_new = P_p_new[d]
            P_w_r_new = P_w_new[r]
            P_p_r_new = P_p_new[r]

            R_d_new = []
            for i in range(len(d)):
                R = self.R_sc(
                    self.p_error, cnr_w_d[i]*P_w_d_new[i], self.m, cnr_p_d[i]*P_p_d_new[i], self.s)
                R_d_new.append(R)

            R_r_new = np.log2(1+factor_w_r*P_w_r_new) + \
                np.log2(1+factor_p_r*P_p_r_new)

            R_new = (np.sum(R_d_new)+np.sum(R_r_new))*sub_band
            return 1,P_w, P_p, R_new

        elif d.size == 0:
            return 0,P_w, P_p, R_reuse
        else:
            return 2,self.all_div_PA(cnr_w, tx_power_w, cnr_p, tx_power_p)

    def all_div_PA(self, cnr_w, tx_power_w, cnr_p, tx_power_p):
        N = self.sm.num_subcarriers
        sub_band = self.sm.spacing
        P_w_all = tx_power_w
        P_p_all = tx_power_p

        k_p = np.exp(self.s * norm.ppf(self.p_error) - self.s ** 2 / 2)
        k_w = gamma.ppf(self.p_error, self.m, scale=1 / self.m)
        #print('k:',k_w, k_p)
        factor_w = cnr_w * k_w
        factor_p = cnr_p * k_p

        p_w_0 = np.ones(N)*P_w_all/N
        p_p_0 = np.ones(N)*P_p_all/N

        snr_w0 = p_w_0*factor_w
        snr_p0 = p_p_0*factor_p
        R_W0 = np.log2(1+snr_w0)
        R_P0 = np.log2(1+snr_p0)
        snr_div0 = np.power(2, R_W0+R_P0)-1

        factor_w_ref = 1e3*k_w * N / P_w_all
        factor_p_ref = 1e3*k_p * N / P_p_all

        def div_equation(x):
            X_w = x[0:N]
            X_p = x[N:2*N]
            P_w = np.exp(X_w)
            P_p = np.exp(X_p)
            X_fi = x[2*N:3*N]
            
            snr_div = np.exp(X_fi)
            snr_w = P_w*cnr_w
            snr_p = P_p*cnr_p

            fw = self.p_pdf_hw_2(snr_div/snr_w)
            fp = self.p_pdf_hp_2(snr_div/snr_p)
            Fw = self.p_cdf_hw_2(snr_div/snr_w)
            Fp = self.p_cdf_hp_2(snr_div/snr_p)

            delta_w = np.sum(P_w)-P_w_all-1/factor_w_ref
            delta_p = np.sum(P_p)-P_p_all-1/factor_p_ref
            yw_part = (1+snr_div)*(1+snr_w * fp*Fw/(snr_p*fw*Fp))
            yw = P_w*yw_part/snr_div+delta_w+0.1*X_w+0.5*X_fi
            yp_part = (1+snr_div)*(1+snr_p * fw*Fp/(snr_w*fp*Fw))
            yp = P_p*yp_part/snr_div+delta_p+0.1*X_p+0.5*X_fi

            log_pe = np.log(self.p_error)-np.log(1e-12+Fp)-np.log(1e-12+Fw)
            result = np.concatenate((yw,  yp, log_pe), axis=0)
            return result

        x_w = np.log(p_w_0)
        x_p = np.log(p_p_0)
        x_fi = np.log(snr_div0)
        x0 = np.concatenate(
            (x_w, x_p, x_fi), axis=0)
        result = np.round(fsolve(div_equation, x0), 2)

        P_w = np.exp(result[0:N])
        P_p = np.exp(result[N:2*N])

        P_w = np.clip(P_w, a_min=1e-6, a_max=None)
        P_p = np.clip(P_p, a_min=1e-6, a_max=None)
        P_w = P_w * P_w_all / np.sum(P_w)
        P_p = P_p * P_p_all / np.sum(P_p)

        R_div = []
        for i in range(N):
            R = self.R_sc(
                self.p_error, cnr_w[i]*P_w[i], self.m, cnr_p[i]*P_p[i], self.s)
            R_div.append(R)

        R_div = np.array(R_div)*sub_band
        return P_w, P_p, R_div

    def reuse_PA(self, cnr_w, tx_power_w, cnr_p, tx_power_p):
        N = self.sm.num_subcarriers
        sub_band = self.sm.spacing
        P_w_all = tx_power_w
        P_p_all = tx_power_p

        k_p = np.exp(self.s * norm.ppf(self.p_error) - self.s ** 2 / 2)
        k_w = gamma.ppf(self.p_error, self.m, scale=1 / self.m)
        factor_w = cnr_w * k_w
        factor_p = cnr_p * k_p
        P_w = np.clip((P_w_all + np.sum(1 / factor_w)) /
                      N - 1 / factor_w, a_min=1e-6, a_max=None)
        P_p = np.clip((P_p_all + np.sum(1 / factor_p)) / N -
                      1 / factor_p, a_min=1e-6, a_max=None)
        P_w = P_w * P_w_all / np.sum(P_w)
        P_p = P_p * P_p_all / np.sum(P_p)
        R_w = np.log2(1+factor_w*P_w)
        R_p = np.log2(1+factor_p*P_p)
        R_w_all = np.sum(R_w)
        R_p_all = np.sum(R_p)
        R_reuse = (R_w_all+R_p_all)*sub_band

        return P_w, P_p, R_reuse
    
    def div_PA(self, cnr_w, tx_power_w, cnr_p, tx_power_p):
        flag,P_w,P_p,R = self.adap_PA(cnr_w, tx_power_w, cnr_p, tx_power_p)
        if flag ==2:
            return P_w,P_p,R
        else:
            N = self.sm.num_subcarriers
            sub_band = self.sm.spacing
            R_div = []
            for i in range(N):
                R = self.R_sc(
                    self.p_error, cnr_w[i]*P_w[i], self.m, cnr_p[i]*P_p[i], self.s)
                R_div.append(R)

            R_div = np.sum(np.array(R_div))*sub_band
            return P_w,P_p,R_div
