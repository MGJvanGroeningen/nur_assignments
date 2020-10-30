#!/usr/bin/env python3
#2d

from ex2_module import np, plt, histogram, divide_in_arrays, P, selection_sort_f

with open('galaxies.txt') as f:
    galaxies = np.array([line[:-1] for line in f.readlines()], dtype = np.float64)

#Find the bin with the largest number of galaxies
bins = np.logspace(-4, np.log10(5), 20)
counts = histogram(galaxies, bins)
ind_max = list(counts).index(max(counts))
big_bin_interval = [bins[ind_max], bins[ind_max + 1]]

#Sort the galaxies in the bin containing most galaxies
big_bin_vals = [gal for gal in galaxies if gal >= big_bin_interval[0] and gal <= big_bin_interval[1]]
big_bin_size = len(big_bin_vals)
sort_big_bin_vals = selection_sort_f(big_bin_vals)

#Save the median, 16th and 84th percentile of the galaxies in the bin containing most galaxies
np.savetxt('nur_a2_2d.txt', (sort_big_bin_vals[int(0.16*big_bin_size)], 
                             sort_big_bin_vals[int(0.50*big_bin_size)], 
                             sort_big_bin_vals[int(0.84*big_bin_size)]))

#Randomly divide the array of galaxies in sub arrays, each containing 100 galaxies
haloes = divide_in_arrays(galaxies, 100)

#Find the number of galaxies contained in the bin containing most galaxies for each halo
gals_in_big_bin = [len([gal for gal in halo if gal >= big_bin_interval[0] and gal <= big_bin_interval[1]]) for halo in haloes]

#Create the poisson distribution
xs = np.arange(0,100)
poisson_dist = np.array([P(sum(gals_in_big_bin)/len(gals_in_big_bin), x) for x in xs])

#Plot the results
fig, ax = plt.subplots()
fig.set_size_inches(10, 6)

ax.hist(gals_in_big_bin, bins = np.arange(100))
ax.plot(xs, poisson_dist*100, label = 'poisson')

ax.set_xlabel('Galaxies', fontsize = 13)
ax.set_ylabel('Occurence', fontsize = 13)
ax.set_title('Probability to find a number galaxies (from a 100 galaxy sample) in the largest bin', fontsize = 15)
ax.legend()

plt.savefig('nur_a2_2d.pdf')

