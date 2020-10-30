#!/usr/bin/env python
#b

from ex1_module import np, plt, my_rand, cdf_hernquist, cdf_hernquist_inv, fraction

start = 1
end = 30000
particles = 1000000

#Create a herquist particle distribution from applying the inverse CDF to an array of uniform values
uni = my_rand(0, 1, particles)
distribution = cdf_hernquist_inv(uni)

#Calculate the fraction of enclosed particles for a number of distances
size = 100
fracs = np.zeros((2, size))
distance = np.logspace(np.log10(start), np.log10(end), size)
for ind, frac_border in enumerate(distance):
    #particle distribution
    fracs[0, ind] = fraction(distribution, frac_border)
    #analytic distribution
    fracs[1, ind] = cdf_hernquist(frac_border)
    
fig, ax = plt.subplots()
fig.set_size_inches(10, 6)

ax.plot(distance, fracs[0], label = 'particle distribution')
ax.plot(distance, fracs[1], '--', label = 'analytic distribution')

ax.set_xscale('log')
ax.set_xlabel(r'$r$ (kpc)', fontsize = 13)
ax.set_ylabel('fraction', fontsize = 13)
ax.set_title('Enclosed fraction of particles of a Hernquist profile', fontsize = 15)
ax.legend()

plt.savefig('nur_a2_1b.pdf')
