#!/usr/bin/env python3
"""2a"""

from ex2_module import *

#Create the particles
np.random.seed(121)
positions = np.random.uniform(low=0, high=16, size = (3,1024))

#Calculate the density at the grid points
grid = np.zeros((16, 16, 16))

for particle in positions.T:
    #Find the indices of the grid points that enclose the particle
    inds = box_nns(np.floor(particle + 0.5))
    for ind in inds:
        #Calculate the weight for each grid point
        delta = weight(particle, ind)
        
        #Apply periodicity
        ind[ind == 16] = 0
        ind[ind == -1] = 15
        
        #Add weights to the grid points
        grid[ind[0], ind[1], ind[2]] += delta

#Calculate the density contrast
mean_dens = 1024/16**3
grid = (grid - mean_dens) / mean_dens

zs = [4, 9, 11, 14]

#Plot the density contrast (delta)
fig, axs = plt.subplots(2,2)
fig.set_size_inches(16, 16)
fig.suptitle(r'Density $\delta$', fontsize=20)

for col in range(2):
    for row in range(2):
        ax = axs[row, col]
        plot = ax.imshow(grid[:,:,zs[col + 2*row]])
        ax.set_title('z = {}'.format(zs[col + 2*row]), fontsize = 18)
        fig.colorbar(plot, ax=ax)

plt.savefig('density.pdf')




"""2b"""

H = np.array(grid, dtype='complex')

#Calculate the Fourier transform of the density contrast
for i in range(16):
    for j in range(16):
        H[i,j,:] = Cooley_Tukey_FFT(H[i,j,:])
for j in range(16):
    for k in range(16):
        H[:,j,k] = Cooley_Tukey_FFT(H[:,j,k])
for i in range(16):
    for k in range(16):
        H[i,:,k] = Cooley_Tukey_FFT(H[i,:,k])

#Calculate the Fourier transform of the gravitational potential (in arbitrary units), by dividing by k^2
for i in range(16):
    for j in range(16):
        for k in range(16):
            if i != 0 or j != 0 or k != 0:
                H[i,j,k] = H[i,j,k]/(i**2 + j**2 + k**2) 

#Slices to plot
zs = [4, 9, 11, 14]                

#Plot the logarithm of the Fourier transform of the gravitational potential
fig, axs = plt.subplots(2,2)
fig.set_size_inches(16, 16)
fig.suptitle(r'Logarithm of the Fourier transform of the gravitational potential $\log_{10}(|\tilde{\Phi}|)$', fontsize=20)

for col in range(2):
    for row in range(2):
        ax = axs[row, col]
        plot = ax.imshow(np.log10(np.abs(H[:,:,zs[col + 2*row]])))
        ax.set_title('z = {}'.format(zs[col + 2*row]), fontsize = 18)
        fig.colorbar(plot, ax=ax)
        
plt.savefig('log_ft_of_phi.pdf')

#Calculate the gravitational potential with inverse Fourier transformation
for i in range(16):
    for j in range(16):
        H[i,j,:] = Cooley_Tukey_FFT(H[i,j,:], inverse=True)
for j in range(16):
    for k in range(16):
        H[:,j,k] = Cooley_Tukey_FFT(H[:,j,k], inverse=True)
for i in range(16):
    for k in range(16):
        H[i,:,k] = Cooley_Tukey_FFT(H[i,:,k], inverse=True)

#Plot the gravitational potential
fig, axs = plt.subplots(2,2)
fig.set_size_inches(16, 16)
fig.suptitle(r'Gravitational potential $\Phi$', fontsize=20)

for col in range(2):
    for row in range(2):
        ax = axs[row, col]
        plot = ax.imshow(H[:,:,zs[col + 2*row]].real)
        ax.set_title('z = {}'.format(zs[col + 2*row]), fontsize = 18)
        fig.colorbar(plot, ax=ax)

plt.savefig('phi.pdf')

