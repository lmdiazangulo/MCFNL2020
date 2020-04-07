import numpy as np 
import matplotlib.pyplot as plt 

def gauss(x,sigma,t=0,c=1,A=1):
    return A*np.exp(-((x-c*t)/(sigma*np.sqrt(2)))**2)

x = np.linspace(-200,200,int(2000))

plt.plot(x,gauss(x,10,50))
plt.show()
