#!/usr/bin/env python
#d

from ex1_module import np, hernquist, hernquist_der, ridders, a 

#Parameters for ridders method
h = 7
d = 2
m = 5

#Calculate the derivative with ridder's method and analytically and their difference
approx = ridders(hernquist, 1.2*a, h, d, m)
analytic = hernquist_der(1.2*a)
diff = np.abs(approx - analytic)

np.savetxt('nur_a2_1d.txt', (approx, analytic, diff))
    
