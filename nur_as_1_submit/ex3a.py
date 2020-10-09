#!/usr/bin/env python3

#3a

from ex3_module import *

x = np.arange(1e-4,5,1e-5)

#The Jacobian for spherical integration is r^2 sin(theta),
#so here integration over every angle gives 4 * pi * x^2
A = N_sat/simpson(4 * np.pi * x**2 * n(x), 1e-5)
np.savetxt('3a_output.txt', [A], fmt ='%1.10f')
