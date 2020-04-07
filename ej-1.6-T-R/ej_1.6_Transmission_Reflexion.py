# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 13:56:49 2020

@author: María Pedrosa Bustos (mpedrosab@gmail.com)
"""

# =============================================================================
# Exercise 1.6: T and R of different panels
# Compute and plot T(!) and R(!) for the following panels
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import speed_of_light, epsilon_0, mu_0




fIni = 1e3
fEnd = 1000e9
N = 1e2 + 1
  
epsilon_r=0.5
sigma=10e3
d=10e-3
mu=1.0*mu_0

def _denominator(phi):         #Computes denominator. _ implies it is private in c++. In python does not implies anything, just for people to read it
    return phi[0][0]*eta_0+phi[0][1]+phi[1][0]*eta_0**2+phi[1][1]*eta_0

omega= np.linspace(fIni,fEnd,N) *2*np.pi        #from freq to wavelength
eta_0 = np.sqrt(mu_0/epsilon_0)
epsilon_c= epsilon_r*epsilon_0+sigma /(np.complex(0,1)*omega)
eta=np.sqrt(mu/epsilon_c)
gamma=np.complex(0,1)*omega*np.sqrt(mu*epsilon_c)
gd=gamma*d

#Matriz transmision
phi=np.array([[np.cosh(gd), eta*np.sinh(gd)], \
               [np.sinh(gd)/eta, np.cosh(gd)]])

#Coeficiente transmisión
T = 2*eta_0/(_denominator(phi))
R = (phi[0][0]*eta_0+phi[0][1]-phi[1][0]*eta_0**2-phi[1][1]*eta_0)/(_denominator(phi))

plt.plot(omega/2/np.pi,np.abs(T),label="T")
plt.plot(omega/2/np.pi,np.abs(R),label="R")
plt.title('$\epsilon_{r}$=%.3f; $\sigma$=%.3e S/m; d=%.3e m' %(epsilon_r,sigma,d))
plt.legend()
