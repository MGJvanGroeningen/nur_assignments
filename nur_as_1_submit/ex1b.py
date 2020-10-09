#!/usr/bin/env python
#1b

from ex1_module import *

Z = 0.5 #new metallicity
n_H = 1e-4 #hydrogen density
tot_Lam = tot_coolrate(metal_free_coolrates, 
                       metal_coolrates, 
                       Electron_densities_over_n_h, 
                       Electron_densities_over_n_h_sol, 
                       Z) #new total coolrate Lambda

#Plots the normalized net coolring rate as function of T at 50% solar metallicity and density 0.0001 cm^-3
#for different values of z. Saves the 101 plots to directory plots1b
for z in np.around(np.linspace(0, 8.989, 101), decimals = 4):
    fig, ax = plt.subplots()
    fig.set_size_inches(10,6)
    ax.plot(Ts, [interp_coolrate(z, T, n_H, tot_Lam) for T in Ts])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(1e-24, 1e-21)
    ax.set_xlim(1e4, 1e9)
    ax.set_title('Normalized net cooling rate at {:.0%}'.format(Z) + r' solar metallicity, $n_{H}$ = 0.001 $\mathrm{cm^{-3}}$' + r' and redshift $z$ = {:.4f} '.format(z))
    ax.set_xlabel(r'$T$ (K)')
    ax.set_ylabel(r'$| \Lambda \; / \; n_{H}^{2} | \; \mathrm{(erg \; cm^{-3} \; s^{-1})}$')
    plt.savefig('./plots1b/snapz_{:.4f}.png'.format(z))
    plt.close(fig)
    
    
