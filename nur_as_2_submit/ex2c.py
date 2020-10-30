#!/usr/bin/env python3
#2c

from ex2_module import np, plt, rand_pick, within_r, selection_sort_f

with open('galaxies.txt') as f:
    galaxies = np.array([line[:-1] for line in f.readlines()], dtype = np.float64)

#Pick 100 random galaxies from the list with 10000 galaxies
gal_picks = rand_pick(galaxies, 100)

#Sort the 100 galaxies
sorted_gal_picks = selection_sort_f(gal_picks)

#Calculate the number of galaxies within a certain radius for an logarithmic array of radii
rs = np.logspace(-4, np.log10(5), 100)
counts = [within_r(sorted_gal_picks, r) for r in rs]

#Plot the results
fig, ax = plt.subplots()
fig.set_size_inches(10, 6)

ax.plot(rs, counts)
ax.set_xscale('log')
ax.set_xlabel(r'$r$ (kpc)', fontsize = 13)
ax.set_ylabel('galaxies', fontsize = 13)
ax.set_title(r'Number of galaxies within $r$', fontsize = 15)

plt.savefig('nur_a2_2c.pdf')

