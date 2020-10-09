#!/usr/bin/env python
import h5py
import numpy as np
import glob
import matplotlib.pyplot as plt

#Fetching the files in CoolingTables
cooltab_loc = './CoolingTables/'
filenames = sorted([file[16:28] for file in glob.glob(cooltab_loc + 'z_*.*.hdf5')])[:-1]
files = [h5py.File(cooltab_loc + filename, 'r') for filename in filenames]

#From the coolingtable files, store the variables used for interpolating
#zs contains all redshift values, Ts contains the temperature bins,
#n_Hs contains the hydrogen density bins and metals contains the metal names
zs = np.sort(np.array([float(file[18:23]) for file in glob.glob(cooltab_loc + 'z_*.*.hdf5')]))
Ts = np.array(files[0]['Metal_free/Temperature_bins'])
n_Hs = np.array(files[0]['Metal_free/Hydrogen_density_bins'])
metals = [name.decode('UTF-8') for name in files[0]['Header/Species_names']]

#Create the spaces in which we want to interpolate, for example H_He_coolrates contains
#the metal free net cooling rates for different z, T and n_H and has shape (49, 352, 81)
metal_free_coolrates = np.stack([np.array(f['Metal_free/Net_Cooling'])[2] for f in files])
Electron_densities_over_n_h = np.stack([np.array(f['Metal_free/Electron_density_over_n_h'])[2] for f in files])
Electron_densities_over_n_h_sol = np.stack([np.array(f['Solar/Electron_density_over_n_h']) for f in files])
metal_coolrates = np.stack([np.stack([np.array(f[metal + '/Net_Cooling']) for f in files]) for metal in metals])

#Caluculates the total normalized net cooling rate
def tot_coolrate(H_He_Lam, metal_Lams, ne_nH, ne_nH_sol, Z):
    return H_He_Lam + sum(metal_Lams) * (ne_nH/ne_nH_sol) * Z
    
def interp_3D(x, x0, y0):
    #adj will contain the indices of the x's 6 adjacent points in x0 in 3D space
    adj = np.zeros((3,2), dtype=np.int16)
    for d in range(3):
        for ind, i in enumerate(x0[d]):
            if i >= x[d]:
                adj[d,0] = ind - 1
                adj[d,1] = ind
                break
    #x_d will express the (normalized) distance between x and and its adjacent points 
    x_d = np.zeros(3)
    for d in range(3):
        x_d[d] = (x[d] - x0[d][adj[d,0]])/(x0[d][adj[d,1]] - x0[d][adj[d,0]])
    #cxx will have the same value as in x in one dimension
    cxx = np.zeros(4)
    for i, j in [[0,0],[0,1],[1,0],[1,1]]:
        cxx[2*i + j] = y0[adj[0,0], adj[1,i], adj[2,j]]*(1 - x_d[0]) + y0[adj[0,1], adj[1,i], adj[2,j]] * x_d[0]
    #cx will have the same value as in x in two dimensions
    cx = np.zeros(2)
    cx[0] = cxx[0] * (1 - x_d[1]) + cxx[2] * x_d[1]
    cx[1] = cxx[1] * (1 - x_d[1]) + cxx[3] * x_d[1]
    #c will contain the interpolated value in x  
    c = cx[0] * (1 - x_d[2]) + cx[1] * x_d[2]
    return c

#Interpolates the cooling rate
def interp_coolrate(z, T, n_H, coolrates):    
    return interp_3D([z, T, n_H], [zs, Ts, n_Hs], coolrates)
    
    
