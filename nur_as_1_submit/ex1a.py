#!/usr/bin/env python
#1a

from ex1_module import *

z = 3 #redshift
Z = 0.25 #metallicity
n_H = [1., 1e-2, 1e-4, 1e-6] #hydrogen density
tot_Lam = tot_coolrate(metal_free_coolrates,
                       metal_coolrates,
                       Electron_densities_over_n_h, 
                       Electron_densities_over_n_h_sol,  
                       Z) #total coolrate Lambda

#Plots the normalized net coolring rate at redshift z = 3, 25% solar metallicity and a few values for the hydrogen density
fig, ax = plt.subplots()
fig.set_size_inches(10,6)
for i in range(4):
    ax.plot(Ts, [interp_coolrate(z, T, n_H[i], tot_Lam) for T in Ts], 
            label = r'$n_{H}$' + ' = {}'.format(n_H[i]) + r' $\mathrm{cm^{-3}}$')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_title(r'Normalized net coooling rate at redshift $z$ = {} and {}% solar metallicity'.format(z, int(100*Z)))
ax.set_xlabel(r'$T$ (K)')
ax.set_ylabel(r'$| \Lambda \; / \; n_{H}^{2} | \; \mathrm{(erg \; cm^{-3} \; s^{-1})}$')
ax.legend()
ax.set_ylim(1e-24, 1e-21)
ax.set_xlim(1e4, 1e9)
plt.savefig('1a_plot.pdf')
    
