# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 13:48:21 2020

@author: María Pedrosa Bustos (mpedrosab@gmail.com)
"""

# =============================================================================
# 
# Exercise 1.2: Fast Fourier transform of a Gaussian pulse
# Implement a python program which, using the module numpy.fft computes
# the fast Fourier transform (FFT), F(!), of a Gaussian
# =============================================================================
#%%
import numpy as np
import matplotlib.pyplot as plt

N=1e3+1  #Num puntos. El +1 es para que el linspace de numeros bonitos
t=np.linspace(0,500e-6,N)
s0=10e-6
t0=10*s0

#Mi función que quiero transformar
f=np.exp(-np.power((t-t0),2)/(2.*np.power(s0,2)))

#Transformada de Fourier
Fourier=np.fft.fft(f)

#Calculo la frecuencia
frec=np.fft.fftfreq(int(N))

#La frecuencia no  la da normalizado a la frecuencia, lo tengo que hacer a mano
frec=frec/(t[1]-t[0])

#Las frecuencias no las ordena bien, entoces sale una línea debajo. Lo reordeno
frec=np.fft.fftshift(frec)
Fourier=np.fft.fftshift(Fourier)

#Ploteo
plt.plot(frec,np.abs(Fourier))