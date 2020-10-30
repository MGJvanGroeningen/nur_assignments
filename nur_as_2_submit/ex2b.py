#!/usr/bin/env python3
#2b

from ex2_module import np, plt, reject_sample, p, N, histogram

#Retrieve the minimum obtained in 2a
with open('nur_a2_2a.txt', 'r') as f:
    x_max = np.array([(f.readline()[:-1])]).astype(np.float64)[0]

start = 0.0001
end = 5
size = 10000
N_sat = 100

#Generate points distributed according to the satellite probability distribution function with rejection sampling
galaxies = reject_sample(p, start, end, size, norm = N_sat/N(x_max))
np.savetxt('galaxies.txt', galaxies)

#Create the bins
bins = np.logspace(-4,np.log10(5), 20)
counts = np.array(histogram(galaxies, bins), dtype = float)

#Weight the bins by dividing by their width
#and divide by the total number of points to normalize the probability distribution
for i in range(len(counts)):
    counts[i] = counts[i]/((bins[i+1] - bins[i])*size)

#Plot the results
fig, ax = plt.subplots()
fig.set_size_inches(10, 6)

#plot the bins
ax.hist(bins[:-1], bins = bins, weights = counts)

#plot the satellite distribution function
x = np.logspace(np.log10(start), np.log10(5), 1000)
ax.plot(x, p(x))

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(1e-4, 10)
ax.set_xlabel(r'$r$ (kpc)', fontsize = 13)
ax.set_ylabel(r'$p(x)$', fontsize = 13)
ax.set_title('Probability distribution of satellite galaxies', fontsize = 15)

plt.savefig('nur_a2_2b.pdf')

