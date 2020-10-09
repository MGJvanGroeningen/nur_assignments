#!/usr/bin/env python3

#3c

from ex3_module import np, P

#Inputs for the Poisson distribution
lams = np.array([1, 5, 3, 2.5, 101], dtype = np.float64)
ks = np.array([0, 10, 21, 40, 200], dtype = np.int64)

np.savetxt('3c_output.txt', [P(lam,k) for lam, k in zip(lams, ks)], fmt ='%1.5e' )
