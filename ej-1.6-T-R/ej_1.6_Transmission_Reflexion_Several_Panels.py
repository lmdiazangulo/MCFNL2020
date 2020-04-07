# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 13:56:49 2020

@author: María Pedrosa Bustos (mpedrosab@gmail.com)
"""

# =============================================================================
# Exercise 1.6: T and R of different panels
# Class for a panel
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import speed_of_light, epsilon_0, mu_0





class Panels():
    
    def __init__():
        self.fIni = 1e3
        self.fEnd = 1000e9
        self.N = 1e2 + 1
      
    self.epsilon_r=epsilon_r
    self.sigma=10e3
    self.d=10e-3
    mu=1.0*mu_0
    
    def _denominator(self, phi):         #Computes denominator. _ implies it is private in c++. In python does not implies anything, just for people to read it
        return phi[0][0]*eta_0+phi[0][1]+phi[1][0]*eta_0**2+phi[1][1]*eta_0
    
    
    self.omega= np.linspace(fIni,fEnd,N) *2*np.pi        #from freq to wavelength
    self.eta_0 = np.sqrt(mu_0/epsilon_0)
    self.epsilon_c= epsilon_r*epsilon_0+sigma /(np.complex(0,1)*omega)
    self.eta=np.sqrt(mu/epsilon_c)
    self.gamma=np.complex(0,1)*omega*np.sqrt(mu*epsilon_c)
    self.gd=gamma*d
    
    #Matriz transmision
    phi=np.array([[np.cosh(gd), eta*np.sinh(gd)], \
                   [np.sinh(gd)/eta, np.cosh(gd)]])
    
    #Coeficiente transmisión
    def T(self, omega)
    T = 2*eta_0/(_denominator(phi))
    R = (phi[0][0]*eta_0+phi[0][1]-phi[1][0]*eta_0**2-phi[1][1]*eta_0)/(_denominator(phi))

plt.plot(omega/2/np.pi,np.abs(T),label="T")
plt.plot(omega/2/np.pi,np.abs(R),label="R")
plt.title('$\epsilon_{r}$=%.3f; $\sigma$=%.3e S/m; d=%.3e m' %(epsilon_r,sigma,d))
plt.legend()
