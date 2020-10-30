#!/usr/bin/env python
#e

from ex1_module import np, hernquist, hernquist_der, newton, cdf_hernquist, M_dm 

crit_dens = 150
delta = [200, 500]

#Define functions turns this problem in a root finding problem
def func_200(r):
    return hernquist(r) - delta[0] * crit_dens

def func_500(r):
    return hernquist(r) - delta[1] * crit_dens

R = np.zeros(2)
M = np.zeros(2)

start = 0.1
iterations = 100
max_err = 1e-10

#Calculate the root with newton raphson method to find the R's
R[0] = newton(func_200, hernquist_der, start, max_err, iterations)[0]
R[1] = newton(func_500, hernquist_der, start, max_err, iterations)[0]

#Multiply the mass fraction with the total mass to get the M's
M = cdf_hernquist(R) * M_dm

np.savetxt('nur_a2_1e1.txt', R, fmt = '%f')
np.savetxt('nur_a2_1e2.txt', M)
    
