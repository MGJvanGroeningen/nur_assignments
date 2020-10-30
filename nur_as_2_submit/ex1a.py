#!/usr/bin/env python
#a

from ex1_module import np, plt, my_rand, pearson

#Initialization of the random variable
#Every time RNG function 'my_rand' runs, it retrieves the random variable from a text file and sets it when its doen
#This text file containing the random variable 'rand.txt' is first created here. 
rand = 1
print('Start with seed = {}'.format(rand))
np.savetxt('rand.txt', np.array([rand]), fmt = '%i')

#Create three different sized arrays containing random numbers between 0 and 1
rand_arrs = [my_rand(size = size) for size in [1000, 1000000, 100000]]

#Calculate the Pearson correlation coefficient for nearest and next nearest elements
r_xi_xi1 = pearson(rand_arrs[2][:-1], rand_arrs[2][1:])
r_xi_xi2 = pearson(rand_arrs[2][:-2], rand_arrs[2][2:])

#Save the Pearson coefficients
np.savetxt('nur_a2_1a.txt', (r_xi_xi1, r_xi_xi2), fmt = '%1.6f')

fig, ax = plt.subplots()
fig.set_size_inches(10, 6)

#Sequential random numbers plotted against each other
ax.scatter(rand_arrs[0][:-1], rand_arrs[0][1:])

ax.set_xlabel(r'$x_{i}$', fontsize = 13)
ax.set_ylabel(r'$x_{i + 1}$', fontsize = 13)
ax.set_title('Sequential random numbers plotted against each other', fontsize = 15)

plt.savefig('nur_a2_1a1.pdf')

fig, ax = plt.subplots()
fig.set_size_inches(10, 6)

#Histogram of the distribution of random numbers between 0 and 1
ax.hist(rand_arrs[1], 20)
ax.hlines((50000 + np.sqrt(50000), 50000 - np.sqrt(50000)), 0, 1, color = 'black', label = r'$\mu \pm \sigma$')

ax.set_ylim(49000, 51000)
ax.set_xlabel(r'$x_{i}$', fontsize = 13)
ax.set_ylabel('N', fontsize = 13)
ax.set_title('Distribution of random numbers from RNG function', fontsize = 15)
ax.legend()

plt.savefig('nur_a2_1a2.pdf')
