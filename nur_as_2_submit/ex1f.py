#!/usr/bin/env python
#f

from ex1_module import np, plt, hernquist_pot, downhill_simplex

size = 150
distances = np.zeros(size)
iterations = np.arange(size)

#Find the minimum of the potential with the downhill simplex method
for i in iterations:
    #Calculate the new best point
    best_point = downhill_simplex(hernquist_pot, [-1000, -200], 0, i)
    #Save the distance from the minimum
    distances[i] = np.sqrt((best_point[0] - 1.3)**2 + (best_point[1] - 4.2)**2)

fig, ax = plt.subplots()
fig.set_size_inches(10, 6)

ax.plot(iterations, distances)
ax.set_yscale('log')
ax.set_xlabel('iterations', fontsize = 13)
ax.set_ylabel('distance (kpc)', fontsize = 13)
ax.set_title('Distance to the minimum of the 2D Hernquist potential', fontsize = 15)

plt.savefig('nur_a2_1f.pdf')
    
