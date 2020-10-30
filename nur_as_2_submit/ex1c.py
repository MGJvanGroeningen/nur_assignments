#!/usr/bin/env python
#c

from ex1_module import np, plt, my_rand, cdf_hernquist_inv, mplot3d

#Create a list of r, theta and phi values of particles from a Hernquist distribution
r = cdf_hernquist_inv(my_rand(0,1,1000))
theta = my_rand(0, np.pi, 1000)
phi = my_rand(0, 2 * np.pi, 1000)

#Transform to Cartesian coordinates
x = r * np.sin(theta) * np.cos(phi)
y = r * np.sin(theta) * np.sin(phi)
z = r * np.cos(theta)

#Plot the points in 3D
fig = plt.figure()
fig.set_size_inches(10,6)
ax = plt.axes(projection='3d')

ax.scatter3D(x,y,z)

ax.set_xlabel('x (kpc)', fontsize = 13)
ax.set_ylabel('y (kpc)', fontsize = 13)
ax.set_zlabel('z (kpc)', fontsize = 13)

ax.set_xlim(-2000, 2000)
ax.set_ylim(-2000, 2000)
ax.set_zlim(-2000, 2000)
ax.set_title('3D particle distribution of a Hernquist profile', fontsize = 15)

plt.savefig('nur_a2_1c1.pdf')

#Plot theta against phi
fig, ax = plt.subplots()
fig.set_size_inches(10,6)

ax.scatter(theta, phi)

ax.set_xlabel(r'$\theta$', fontsize = 13)
ax.set_ylabel(r'$\phi$', fontsize = 13)
ax.set_title(r'Correlation between $\theta$ and $\phi$', fontsize = 15)

plt.savefig('nur_a2_1c2.pdf')

    
