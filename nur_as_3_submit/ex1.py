#!/usr/bin/env python

"""This file contains answers for 1a, 1b and 1c"""

print('Running pipeline...')

from ex1_module import *

print('Pipeline complete!')

"""Initialize data products"""

#Initial parameters
par0 = [1.8, 0.8, 2.8]

#Parameters a, b, c for the best fitting model to the data from each file (for each distribution) 
pars = np.zeros((2, 5, 3))

#Gaussian and Poisson statistic measure for each file
chi2s = np.zeros(5)
min_log_liks = np.zeros(5)

#G test statistic and Q values for each file and distribution 
G_tests = np.zeros((2, 5))
Qs = np.zeros((2, 5))


"""Calculate data products"""


dist_name = ['Gaussian', 'Poisson']

if os.path.exists("nur_a3_1_output.txt"): 
    os.remove("nur_a3_1_output.txt")
    
iterations = 100

#File
for m in range(5):
    #Distribution
    print('Calculating best fits for file: {} ...'.format(filenames[m]))
    for dist in range(2):  
        #Define the N_sat belonging to the correct file (used in functions that calculate the model)
        N_sat = N_sats[m]
        n_haloes = n_haloes_list[m]
        
        #Calculate the best parameters and corresponding chi squared for the Gaussian approach
        if dist == 0:
            pars[dist, m], chi2s[m] = levenberg_marquardt(counts[m], par0, 0, iterations, N_sat, n_haloes)
        
        #Calculate the best parameters and corresponding min log likelihood for the Poisson approach
        else:
            #Define new function that only takes the parameters as variable
            def min_ln_lik_poisson_new(par):
                return min_ln_lik_poisson(counts[m], model_vals(par, N_sat, n_haloes))
            pars[dist, m] = downhill_simplex(min_ln_lik_poisson_new, par0, 0, iterations)
        
        #Calculate the model with the best found parameters
        model = model_vals(pars[dist, m], N_sat, n_haloes)
        
        #Statistic for Poisson (the one for Gaussian is given by the Levenberg-Marquardt function)
        if dist == 1:
            min_log_liks[m] = min_ln_lik_poisson(counts[m], model)
        
        #Evaluate model with a Q value from a G-test
        G_tests[dist, m], df = G_test(counts[m], model)
        Qs[dist, m] = Q_value(G_tests[dist, m], df)
        
        print('{} fit done!'.format(dist_name[dist]))
        
        #Also print results to a file
        with open("nur_a3_1_output.txt", "a") as f:
            if dist == 0:
            	print('File: {}'.format(filenames[m]), file=f)
            print(' ', file=f)
            print('Approach: {}'.format(dist_name[dist]), file=f)
            print(' ', file=f)
            print('Parameters:', file=f)
            print('a = {}'.format(pars[dist, m, 0]), file=f)
            print('b = {}'.format(pars[dist, m, 1]), file=f)
            print('c = {}'.format(pars[dist, m, 2]), file=f)
            print(' ', file=f)
            print('Statistics:', file=f)
            print('Degrees of freedom = {}'.format(df), file=f)
            if dist == 0:
                print('Chi squared = {}'.format(chi2s[m]), file=f)
            else:
                print('Min log likelihood = {}'.format(min_log_liks[m]), file=f)
            print('G = {}'.format(G_tests[dist, m]), file=f)
            print('Q = {}'.format(Qs[dist, m]), file=f)
            if dist == 1:
            	print(' ', file=f)
            	print(30*'=', file = f)
            	print(' ', file=f)
            
"""Plot the data products"""            
            
#Radius values of the bins
rs = np.array([(bins[i+1] + bins[i])/2 for i in range(len(bins[:-1]))])

#Approximate masses of the haloes per file (used in plot title)
mass = [r'10^{11}', r'10^{12}', r'10^{13}', r'10^{14}', r'10^{15}']

#Divide the counts by the total number of haloes to get the average number of satellites per bin per halo
weighted_counts = (counts.T/np.array(n_haloes_list)).T

#Plot the models and weighted counts for the distribution of satellite galaxies 
fig, ax = plt.subplots(3, 2)
fig.set_size_inches(16, 18)
fig.delaxes(ax[2,1])
fig.tight_layout(pad =6.)

for m in range(5): 
    ax[int(m/2), m % 2].hist(bins[:-1], bins=bins, weights=weighted_counts[m], label = 'Data', histtype = 'step', linewidth = 2)
    N_sat = N_sats[m]
    for dist in range(2):
        #Scale the model and uncertainty to get the average number of satellites per bin per halo
        model = model_vals(pars[dist, m], N_sat, n_haloes_list[m])/n_haloes_list[m]
        sigma = np.sqrt(model/n_haloes_list[m])
        if dist == 0:
            ax[int(m/2), m % 2].errorbar(rs, model, np.full(len(rs), sigma), label = dist_name[dist] + ' fit', capsize = 3)
        else:
            ax[int(m/2), m % 2].errorbar(rs, model, np.full(len(rs), sigma), label = dist_name[dist] + ' fit', capsize = 3)
    ax[int(m/2), m % 2].set_xscale('log')
    ax[int(m/2), m % 2].set_yscale('log')
    ax[int(m/2), m % 2].set_ylim(min(weighted_counts[m][weighted_counts[m] != 0])/10, 10*max(weighted_counts[m]))
    ax[int(m/2), m % 2].set_xlabel('x', fontsize = 13)
    ax[int(m/2), m % 2].set_ylabel('Average satellite galaxies per halo', fontsize = 13)
    ax[int(m/2), m % 2].set_title(r'Distribution of satellite galaxies in haloes with $M \approx ' + mass[m] + ' \; M_{\odot} \; h^{-1}$', fontsize = 15)
    ax[int(m/2), m % 2].legend(loc = 'upper left')


plt.savefig('nur_a3_1_plot.pdf')
