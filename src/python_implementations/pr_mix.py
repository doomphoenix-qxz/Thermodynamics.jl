import numpy as np
import scipy.optimize as opt
import scipy.integrate as intg
import matplotlib.pyplot as plt
from scipy.misc import derivative


def R():
    return 8.31446


class pengRobinson:
    def __init__(self, Tc, Pc, omega, Cp_func):
        self.Tc = Tc
        self.Pc = Pc
        self.omega = omega
        self.a = 0.45724 * R() ** 2 * Tc ** 2 / Pc
        self.b = 0.07780 * R() * Tc / Pc
        self.k = 0.37464 + 1.54226 * omega - 0.26992 * omega ** 2
        self.Cp = Cp_func

    def alpha(self, T):
        return (1 + self.k * (1 - np.sqrt(T / self.Tc))) ** 2

    def dAdT(self, T):
        Tr = T / self.Tc
        return (np.sqrt(Tr) - 1) * (self.k - 1) * self.k / (self.Tc * np.sqrt(Tr))

    def Z_rho(self, T, rho):
        t1 = 1 / (1 - self.b * rho)
        t2 = self.a * self.alpha(T) / (R() * T * (1 / rho + 2 * self.b - self.b ** 2 * rho))
        return t1 - t2

    def dZdT(self, T_, rho):
        daldT = self.dAdT(T_)
        return -(self.a * daldT - self.alpha(T_) / T_) / (R() * T_ * (1 / rho + 2 * self.b - self.b ** 2 * rho))

    def Hr(self, T, rho):
        def dZdT(T_, rho):
            return derivative(lambda Ts: self.Z_rho(Ts, rho), T_, dx=1e-6)

        integral = intg.quad(lambda rho_: dZdT(T, rho_) / rho_, 0.001, rho, limit=200)[0]
        return R() * T * (self.Z_rho(T, rho) - 1 - T * integral)

    def H(self, P, T):
        return intg.quad(self.Cp, 50, T)[0] + self.Hr(T, 1 / self.volume(P, T))

    def U(self, P, T):
        return self.H(P, T) - P * self.volume(P, T)

    def volume(self, p, T):
        a_ = self.alpha(T)
        A = self.a * a_ * p / T ** 2 / R() ** 2
        B = self.b * p / R() / T
        coeffs = np.array([1, -(1 - B), (A - 2 * B - 3 * B ** 2), -(A * B - B ** 2 - B ** 3)])
        Z = np.roots(coeffs)
        V = float(Z[0]) * R() * T / p
        return V

    def v_2phase(self, p, T):
        a_ = self.alpha(T)
        A = self.a * a_ * p / T ** 2 / R() ** 2
        B = self.b * p / R() / T
        coeffs = np.array([1, -(1 - B), (A - 2 * B - 3 * B ** 2), -(A * B - B ** 2 - B ** 3)])
        Z = np.roots(coeffs)
        V = Z * R() * T / p
        return V
    
    def vv(self, p, T):
        def solveit(Vg):
            return p - self.pressure(Vg, T)

        return float(opt.fsolve(solveit, R() * T / p))

    def vl(self, p, T):
        def solveit(Vg):
            return p - self.pressure(Vg, T)

        return float(opt.fsolve(solveit, 1.1 * self.b))

    def pressure(self, V, T):
        t1 = R() * T / (V - self.b)
        t2 = self.a * self.alpha(T) / (V ** 2 + 2 * self.b * V - self.b ** 2)
        return t1 - t2

    def temperature(self, p, V):
        return opt.fsolve(lambda T_: p - self.pressure(V, T_), 300)

    def dVdT(self, p, T):
        def derthis(T_):
            return self.volume(p, T_)

        return derivative(derthis, T, dx=1e-6)

    def dpdV(self, V, T):
        def derthis(V_):
            return self.pressure(V_, T)

        return derivative(derthis, V, dx=1e-6)

    def find_pstar(self, T, p_guess):
        # Vspinodal_v = opt.fsolve(lambda Vg: self.dpdV(Vg, T), 0.1)
        # pspinodal_v = self.pressure(Vspinodal_v, T)
        # print(self.dpdV(Vspinodal_v, T))
        # Vspinodal_l = opt.fsolve(lambda Vg: self.dpdV(Vg, T), 1.1*self.b)
        # pspinodal_l = self.pressure(Vspinodal_l, T)
        # print(self.dpdV(Vspinodal_l, T))
        # p_guess = (pspinodal_l + pspinodal_v)/2
        def solve_this(pg):
            V3 = self.v_2phase(pg, T)
            vl = V3[0]
            vi = V3[1]
            vv = V3[2]

            # Verify vapor pressure
            def integrator(V, T_):
                return self.pressure(V, T_) - pg

            piece1 = intg.quad(integrator, vl, vi, args=T)[0]
            piece2 = intg.quad(integrator, vi, vv, args=T)[0]
            areadd = piece1 + piece2
            #print(areadd)
            return areadd

        return opt.fsolve(solve_this, p_guess)
        

    def phi(self, p, T, liqflag=False):
        A_ = self.a * p / (R() * T) ** 2
        B_ = self.b * p / (R() * T)
        if liqflag:
            V_ = self.vl(p, T)
        else:
            V_ = self.vv(p, T)
        Z_ = p * V_ / (R() * T)
        t1 = (Z_ - 1) - np.log(Z_ - B_)
        t2 = A_ / (2 * B_ * np.sqrt(2)) * np.log((Z_ + B_ * (1 + np.sqrt(2))) / (Z_ + B_ * (1 - np.sqrt(2))))
        return np.exp(t1 - t2)


class PR_Mixture:
    def __init__(self, prlist, yi, kij):
        self.yi = yi
        self.pr_ = dict(zip(yi, prlist))
        self.kij = dict(zip(yi, kij))
        b_ = 0
        for y_ in yi:
            b_+= y_*self.pr_[y_].b
        self.b = b_

    def amix(self, T):
        ans = 0
        for yi_ in self.yi:
            for yj in self.yi:
                pr1 = self.pr_[yi_]
                pr2 = self.pr_[yj]
                ans += yi_*yj*np.sqrt(pr1.a*pr1.alpha(T)*pr2.a*pr2.alpha(T))*(1 - self.kij[yi_])
        return ans

    def a(self, T):
        return self.amix(T)

    def comp(self, yi_, T):
        ans = 0
        pr1 = self.pr_[yi_]
        for yj in self.yi:
            pr2 = self.pr_[yj]
            ans += yj * np.sqrt(pr1.a * pr1.alpha(T) * pr2.a * pr2.alpha(T))*(1 - self.kij[yi_])
        return 2*ans

    def volume(self, p, T):
        A = self.a(T) * p / T ** 2 / R() ** 2
        B = self.b * p / R() / T
        coeffs = np.array([1, -(1 - B), (A - 2 * B - 3 * B ** 2), -(A * B - B ** 2 - B ** 3)])
        Z = np.roots(coeffs)
        V = float(Z[0]) * R() * T / p
        return V

    def pressure(self, V, T):
        t1 = R() * T / (V - self.b)
        t2 = self.a(T) / (V ** 2 + 2 * self.b * V - self.b ** 2)
        return t1 - t2

    def find_pstar(self, T, p_guess):
        def solve_this(pg):
            V3 = self.v_2phase(pg, T)
            vl = V3[0]
            vi = V3[1]
            vv = V3[2]

            # Verify vapor pressure
            def integrator(V, T_):
                return self.pressure(V, T_) - pg

            piece1 = intg.quad(integrator, vl, vi, args=T)[0]
            piece2 = intg.quad(integrator, vi, vv, args=T)[0]
            areadd = piece1 + piece2
            # print(areadd)
            return areadd

        return opt.fsolve(solve_this, p_guess)

    def vv(self, p, T):
        def solveit(Vg):
            return p - self.pressure(Vg, T)

        return float(opt.fsolve(solveit, R() * T / p))

    def vl(self, p, T):
        def solveit(Vg):
            return p - self.pressure(Vg, T)

        return float(opt.fsolve(solveit, 1.1 * self.b))

    def ln_phi(self, p, T, index=1, liqflag=False):
        i = index-1
        B = self.b * p / R() / T
        yi = self.yi[i]
        bi = self.pr_[yi].b
        Bi = bi * p / R() / T
        if liqflag:
            V_ = self.vl(p, T)
        else:
            V_ = self.vv(p, T)
        Z_ = p * V_ / (T * R())
        sigaj = self.comp(yi, T)
        A = self.a(T) * p / T ** 2 / R() ** 2
        t1 = (Bi / B) * (Z_ - 1) - np.log(Z_ - B)
        t2 = -(A / (2 * np.sqrt(2) * B))
        t4 = (sigaj / self.a(T) - Bi / B)
        t3 = np.log((Z_ + B * (1 + np.sqrt(2))) / (Z_ + B * (1 - np.sqrt(2))))
        return t1 + t2 * t4 * t3



    def phi(self, p, T, index=1, liqflag=False):
        return np.exp(self.ln_phi(p, T, index, liqflag))

    def fi(self, p, T, index=1, liqflag=False):
        return self.yi[index-1] * p * self.phi1(p, T, index, liqflag)


    def f_mix(self, p, T, liqflag=False):
        return self.f_1(p, T, liqflag) + self.f_2(p, T, liqflag)



def __main__():
    ethane = pengRobinson(305.3,48.72e5,0.1,None)
    propane = pengRobinson(369.8,42.48e5,0.152,None)
    butane = pengRobinson(425.1,37.96e5,0.2,None)
    x = [0.125, (2/3), (1 - 0.125 - 2/3)]
    P = 20e5
    T = 342.422
    prlist = [ethane, propane, butane]
    eqmix = PR_Mixture(prlist, x, [0,0,0])
    phi1 = eqmix.phi(P, T, 1)
    phi2 = eqmix.phi(P,T,2)
    phi3 = eqmix.phi(P,T,3)
    a_alph = eqmix.amix(T)
    sigaj = eqmix.comp(x[0], T)
    print("A_alpha mix: {}".format(a_alph))
    print("Comp mix: {}".format(sigaj))
    print("Fugacity coefficients: {} ethane {} propane {} butane".format(phi1, phi2, phi3))


if __name__ == "__main__":
    __main__()