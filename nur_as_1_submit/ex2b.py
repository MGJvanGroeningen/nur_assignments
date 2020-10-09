#!/usr/bin/env python3
#2b

from ex2_module import *
from ex2a import f

#Create a new f via iterative improvement
new_f = f - LU_solve(w_ss, np.dot(w_ss, f) - w_gs)[2]

#Save the output to a text file
np.savetxt('2b_output_f.txt', new_f, fmt ='%1.10f')
np.savetxt('2b_output_sumf.txt', [sum(new_f)], fmt='%1.15f')

