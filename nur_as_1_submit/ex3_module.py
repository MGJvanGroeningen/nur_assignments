#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import Akima1DInterpolator

#Parameters
a = 2.2
b = 0.5
c = 3.1
N_sat = 100

#Gives the number denisty profile for satellite galaxies
def n(x):
    return N_sat * (x/b)**(a - 3) * np.exp(-(x/b)**c)

#Gives the 10-base logarithm of n(x)
def log10_n(x):
    return np.log10(N_sat * (x/b)**(a-3)) + np.log10(np.e)*(-(x/b)**c)

#Integration function via simpson's method
def simpson(y,h):
    return h*(sum(2*y[::2]) + sum(4*y[1::2]) - y[0] - y[-1])/3

#Function for linear interpolation
def lin_interp(x, x_0, y_0):
    y = []
    for i in x:
        for ind, j in enumerate(x_0):
            if i == j:
                y.append(y_0[ind])
            elif i <= x_0[-1] and i > j and i < x_0[ind + 1]:
                y.append(y_0[ind] + (i - x_0[ind]) * (y_0[ind + 1] - y_0[ind]) / (x_0[ind + 1] - x_0[ind]))
                break
    return np.array(y)

#Formula for Poisson distribution
def P(lam, k):
    #To avoid the overflow from dealing with large numbers calculate the Poisson distribution in log space
    #and then convert back to normal space
    #The factorial in log space is simply the summation of a list with the logarithm of all (positive) integers up to k
    return 10**(k*np.log10(lam) + np.log10(np.e) * (-lam) - sum(np.log10(np.arange(1, k+1))))
    
