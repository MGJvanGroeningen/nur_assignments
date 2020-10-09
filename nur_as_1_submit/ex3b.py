#!/usr/bin/env python3

#3b

from ex3_module import *
from ex3a import x, A

#Transfrom x and n(x) to 10-base log space
log_x = np.log10(x)
log_n = A * log10_n(x)

#Create the sample points and transform these to log space as well
sample_x = np.array([1e-4, 1e-2, 1e-1, 1, 5], dtype=float)
log_sample_x = np.log10(sample_x)
log_sample_n = A * log10_n(sample_x)

#Interpolate log(x) using the log sample points (linear and akima) 
lin = lin_interp(log_x, log_sample_x, log_sample_n)
akima = Akima1DInterpolator(log_sample_x, log_sample_n)(log_x)

#Plots the 10-base logarithm of n(x) and a linear and akima interpolation
fig, ax = plt.subplots()
fig.set_size_inches(10,6)
ax.plot(log_x, lin, label = 'linear')
ax.plot(log_x, akima, label = 'akima')
ax.plot(log_x, log_n, label = 'true')
ax.scatter(log_sample_x, log_sample_n)
ax.set_title('Radial distribution of galaxy satellites')
ax.set_ylabel(r'$\log (n(x))$')
ax.set_xlabel(r'$\log (x)$')
ax.legend()
plt.savefig('3b_plot.pdf')
