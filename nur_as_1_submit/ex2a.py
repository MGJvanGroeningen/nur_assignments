#!/usr/bin/env python3
#2a

from ex2_module import *

#Create the lower and upper matrix, solves the matrix equation and returns f
L, U, f = LU_solve(w_ss, w_gs)

#Save the output to text files
np.savetxt('2a_output_L.txt', L, fmt ='%1.3f')
np.savetxt('2a_output_U.txt', U, fmt ='%1.3f')
np.savetxt('2a_output_f.txt', f, fmt ='%1.10f')
np.savetxt('2a_output_sumf.txt', [sum(f)], fmt='%1.15f')

