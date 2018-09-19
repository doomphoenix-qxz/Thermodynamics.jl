# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 09:59:29 2017

@author: richar56
"""

import numpy as np
import scipy.optimize as opt
import scipy.integrate as intg
import matplotlib.pyplot as plt
from scipy.misc import derivative
from pr_mix import R, pengRobinson, PR_Mixture

class NRTL:
    
    def __init__(self, tau12, tau21, alpha):
        # tau assumed to be a function of temperature specified by the user
        self.tau1 = tau12
        self.tau2 = tau21
        self.alpha= alpha
        self.G1 = lambda T: np.exp(-self.alpha*self.tau1(T))
        self.G2 = lambda T: np.exp(-self.alpha*self.tau2(T))
        
    def ln_g1(self, x1, T):
        x2 = 1-x1
        G1 = self.G1(T)
        G2 = self.G2(T)
        t1 = self.tau1(T)
        t2 = self.tau2(T)
        term1 = t2*G2**2 / (x1 + x2*G2)**2
        term2 = t1*G1    / (x2 + x1*G1)**2
        return x2**2 * (term1+term2)
    
    def ln_g2(self, x1, T):
        x2 = 1-x1
        G1 = self.G1(T)
        G2 = self.G2(T)
        t1 = self.tau1(T)
        t2 = self.tau2(T)
        term1 = t2*G2    / (x1 + x2*G2)**2
        term2 = t1*G1**2 / (x2 + x1*G1)**2
        return x1**2 * (term1+term2)
    
    def g1(self, x1, T):
        return np.exp(self.ln_g1(x1, T))
    
    def g2(self, x1, T):
        return np.exp(self.ln_g2(x1, T))
    
    def lle(self, T, cutoff):
        guesses = np.array([0.30, 0.90])
        def solveit(guess):
            x1a, x1b = guess
            x2a = 1-x1a
            x2b = 1-x1b
            a1_a = np.log(x1a) + self.ln_g1(x1a, T)
            a2_a = np.log(x2a) + self.ln_g2(x1a, T)
            a1_b = np.log(x1b) + self.ln_g1(x1b, T)
            a2_b = np.log(x2b) + self.ln_g2(x1b, T)
            err = np.array([a1_a - a1_b, a2_a - a2_b])
            if x1a >= 0.58:
                err += np.abs(x1a-0.58)*100
            if x1b <= 0.68:
                err += np.abs(x1b-0.68)*100
            return err
        ans = opt.fsolve(solveit, guesses)
        return ans
    
    def GE_RT(self, x1, T):
        x2 = 1-x1
        return self.ln_g1(x1, T)*x1 + self.ln_g2(x1, T)*x2
    
    def GE(self, x1, T):
        return self.GE_RT(x1, T) * R() * T
    
    def G_mix(self, x1, T):
        """
        Note that this function gives DeltaG mixture, not the absolute Gibbs
        Energy, which would require also the Gibbs energies of individual
        components.
        """
        x2 = 1-x1
        return R()*T*(x1*np.log(x1) + x2*np.log(x2)) + self.GE(x1, T)
    
    

def tau12_isooctane_perfluoroheptane(T):
    return 1.420 + 155.5/T

def tau21_isooctane_perfluoroheptane(T):
    return 0.785 + 108.0/T

def bubble_calc(p, x1, prs, gamma_model):
    x2 = 1-x1
    pr1 = prs[0]
    pr2 = prs[1]
    def solveit(guess):
        T,y1 = guess
        err = np.zeros(2)
        y2 = 1-y1
        prmix = PR_Mixture(prs, np.array([y1, y2]), np.array([0,0]))
        Ps1 = pr1.find_pstar(T, p)
        Ps2 = pr2.find_pstar(T, p)
        phis1 = pr1.phi(Ps1, T)
        phis2 = pr2.phi(Ps2, T)
        phi_1 = prmix.phi(p, T, 1)
        phi_2 = prmix.phi(p, T, 2)
        g1 = gamma_model.g1(x1, T)
        g2 = gamma_model.g2(x1, T)
        Poynt1 = np.exp(pr1.vl(p, T) * (p - Ps1) / R() / T)
        Poynt2 = np.exp(pr2.vl(p, T) * (p - Ps2) / R() / T)
        err[0] = y1 * phi_1 * p - x1*g1*Ps1*phis1*Poynt1 
        err[1] = y2 * phi_2 * p - x2*g2*Ps2*phis2*Poynt2
        return err
    if x1 > 0.5:
        y1g = x1-0.4
    else:
        y1g = x1+0.25    
    ans = opt.fsolve(solveit, np.array([283.15, y1g]))
    return ans

perf_hept = pengRobinson(475.65, 1610000, 0.542889, None)
isooctane = pengRobinson(543.8, 2570000, 0.303455, None)
        

alpha_ = 0.4

proj_mix = NRTL(tau12_isooctane_perfluoroheptane, tau21_isooctane_perfluoroheptane, alpha_)
x1_space = np.linspace(0, 1, 101)


#print(proj_mix.lle(273.15-50.0, 0.65))
#y1_space = np.zeros(len(x1_space))
#T_space = np.zeros(len(x1_space))
#for i, x1 in enumerate(x1_space):
#    T_space[i], y1_space[i] = bubble_calc(0.05e5, x1, [isooctane, perf_hept], proj_mix) 
    
xs, y1_space, T_space = np.loadtxt("project2_vledata.csv")


Tm_iso = 165.777
dCp_iso = (1.8285e2 - 1.5429e2) # J/mol*K
Tm_pfl = 221.86
dCp_pfl = (377.4 - 337.41)      # J/mol*K

def psi1(T):
    Hm = 9196                     # J/mol
    return np.exp(-Hm * (1/T - 1/Tm_iso)/R() + dCp_iso * (Tm_iso/T - 1 - np.log(Tm_iso/T))/R())

def psi2(T):
    Hm = 6948                     # J/mol
    return np.exp(-Hm * (1/T - 1/Tm_pfl)/R() + dCp_pfl * (Tm_pfl/T - 1 - np.log(Tm_pfl/T))/R())


def SLE(x1, gamma_model):
    def solveit(guess):
        T = guess
        err = np.zeros(1)
        err[0] = np.log(x1) + gamma_model.ln_g1(x1, T) - np.log(psi1(T))
        #err[0] = np.log(1-x1) + gamma_model.ln_g2(x1, T) - psi2(T)
        return err
    ans = opt.fsolve(solveit, [200])
    return ans

def SLE2(x1, gamma_model):
    def solveit(guess):
        T = guess
        err = np.zeros(1)
        #err[0] = np.log(x1) + gamma_model.ln_g1(x1, T) - psi1(T)
        err[0] = np.log(1-x1) + gamma_model.ln_g2(x1, T) - np.log(psi2(T))
        return err
    ans = opt.fsolve(solveit, [200])
    return ans


Ts3 = np.zeros(len(x1_space))
Ts4 = np.zeros(len(x1_space))

for i, x1 in enumerate(x1_space):
    Ts4[i] = SLE(x1, proj_mix)
    Ts3[i] = SLE2(x1, proj_mix)


T_lle = T_space[50]
Ts2 = np.linspace(T_lle, Ts3[40], 101)
xa1_ = np.zeros(len(Ts2))
xa2_ = np.zeros(len(Ts2))

for i, T_ in enumerate(Ts2):
    xa1_[i], xa2_[i] = proj_mix.lle(T_, 0.65)
    
def get_e_pt():
    def solveit(guess):
        x1, T = guess
        err = np.zeros(2)
        err[0] = psi1(T)/proj_mix.g1(x1, T) + psi2(T)/proj_mix.g2(x1, T) - 1
        err[1] = x1*proj_mix.g1(x1, T) - psi1(T)
        return err
    ans = opt.fsolve(solveit, np.array([0.5, 170]))
    
    return ans

ep = get_e_pt()
print(ep)
    
plt.figure(figsize=(8,10))
plt.plot(x1_space, np.ones(len(x1_space))*T_lle, '--m', linewidth=1)
plt.plot(x1_space, T_space)
plt.plot(y1_space, T_space)
plt.plot(x1_space[0:86], np.ones(len(x1_space[0:86])) * Ts3[40], '--', color="xkcd:crimson",linewidth=1)
plt.plot(x1_space[0:42], Ts3[0:42])
plt.plot(x1_space[85:100], Ts3[85:100])

plt.plot(x1_space, np.ones(len(x1_space))*ep[1], '--y', linewidth=1)
plt.plot(x1_space[98:], Ts4[98:], linewidth=2)


plt.plot(xa1_, Ts2)
plt.plot(xa2_, Ts2)
plt.xlim(0,1)
plt.ylim(150, 300)

saveit = np.array([x1_space, y1_space, T_space])

print(proj_mix.g2(0.5,200.0))
print(psi1(200.0))
print(psi2(200.0))
plt.xlabel("Mole fraction of isooctane")
plt.ylabel("Temperature (K)")
plt.title("Phase diagram for isooctane (1) and perfluoroheptane (2)")
plt.savefig("project_diagram.png")

