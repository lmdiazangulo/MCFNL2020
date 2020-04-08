import numpy as np 
import matplotlib.pyplot as plt 

def gaussian(x_array, t, A, c, sigma):
    expon = (x_array - c*t) / sigma / np.sqrt(2)
    return A*np.exp(-np.power(expon, 2))


x = np.linspace(-20.0, 100.0, 1000)
plt.plot(x , gaussian(x, 15, 2, 2, 5))
plt.show()
