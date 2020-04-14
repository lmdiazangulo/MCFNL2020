import numpy as np
import scipy.constants as sp

N = int(1e6)
a = 1e-6

# Product and then sum 
aProdSum = 0
for i in range(N):
    aProdSum += a * sp.mu_0
print ("Prod Sum diff: %f\n", aProdSum - sp.mu_0)

# Sum and then product
aSumProd = 0
for i in range(N):
    aSumProd += a
aSumProd *= sp.mu_0
print ("Sum Prod diff: %f\n", aSumProd - sp.mu_0)

# Rounding
for i in range (20):
    print ( i, i/10, i*.1 )

print("END")