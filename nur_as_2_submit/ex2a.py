#!/usr/bin/env python3
#2a

from ex2_module import np, bracket_finder, golden_search, N_neg, N

#Find a bracket which contains the minimum (0.05 and 0.1 were considered reasonable starting values after plotting N(x))
bracket = bracket_finder(N_neg, 0.05, 0.1)

#Find the minimum in the bracket
x_max = golden_search(N_neg, bracket, 1e-9)[0]

np.savetxt('nur_a2_2a.txt', (x_max, N(x_max)))

