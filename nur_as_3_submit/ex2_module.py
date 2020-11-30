#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

#Nearest neighbour indices
nn_inds = np.array([[0,0,0],
                    [1,0,0],
                    [0,1,0],
                    [0,0,1],
                    [1,1,0],
                    [0,1,1],
                    [1,0,1],
                    [1,1,1]], dtype=np.int8)

#Calculates the weight for a particle and grid point position
def weight(pos1, pos2):
    return ((1 - np.abs(pos1[0] - (pos2[0] + 0.5))) * 
            (1 - np.abs(pos1[1] - (pos2[1] + 0.5))) * 
            (1 - np.abs(pos1[2] - (pos2[2] + 0.5))))

#Finds the indices of the grid points that enclose position pos
def box_nns(pos):
    inds = np.zeros((8, len(pos)), dtype = np.int8)
    for i in range(8):
        inds[i] = np.array(pos, dtype = np.int8) - nn_inds[i]
    return inds

#Discrete Fourier transform routine
def DFT_routine(x, start, N_j, inverse=False):
    N_j2 = int(N_j/2)
    
    #Recursive loops
    if N_j > 2:
        DFT_routine(x, start, N_j2, inverse)
        DFT_routine(x, start + N_j2, N_j2, inverse)
    
    #Add a minus sign when taking the inverse
    if inverse:
        theta = -2 * np.pi / N_j
    else:
        theta = 2 * np.pi / N_j
    
    #Combine values from the array to make the Fourier transform
    alpha = 2 * np.sin(theta / 2)**2
    beta = np.sin(theta)
    for k in range(start, start + N_j2):
        t = x[k]
        phi = (k - 1)*theta
        cos_k_theta = (1 - alpha) * np.cos(phi) - beta * np.sin(phi)
        sin_k_theta = (1 - alpha) * np.sin(phi) + beta * np.cos(phi)
        W_k = cos_k_theta + 1j * sin_k_theta
        x[k] = t + W_k * x[k + N_j2]
        x[k + N_j2] = t - W_k * x[k + N_j2]

#Caculates the (inverse) Fourier transform of array x with the Cooley Tukey algorithm
def Cooley_Tukey_FFT(x, inverse=False):
    x = np.array(x, dtype='complex')
    
    #Padding
    N = len(x)
    k = 0
    while True:
        #Find the nearest power of two that is larger than N
        padding = int((2**k) - N)
        if padding < 0:
            k += 1
        else:
            #Add zeros up to a power of two
            x = np.concatenate((x, np.zeros(padding, dtype='complex')))
            break
    
    #Bit-reversing indices
    N = len(x)
    bits = int(np.log2(N))
    indices = np.array([int('{:0{bits}b}'.format(k, bits=bits)[::-1], 2) for k in range(N)])
    x = np.array([x[ind] for ind in indices])
    
    #Calculate FFT
    DFT_routine(x, 0, N, inverse)
    
    #Divide by N when taking the inverse
    if inverse:
        x = x / N
        
    #Remove the padding
    x = x[:len(x) - padding]
    return x

