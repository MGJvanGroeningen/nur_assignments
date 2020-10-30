#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

G = 4.3009e-6 #kpc M_sol^(-1) (km/s)^2
M_dm = 1e12 # M_sol
a = 80 # kpc
    
#Retrieve the random variable from rand.txt
def get_rand():
    with open('rand.txt', 'r') as f:
        rand = int(f.readline()[:-1])
    return rand

#Write the random variable to rand.txt  
def set_rand(rand):
    np.savetxt('rand.txt', np.array([rand]), fmt = '%i')

#Random number generator
def my_rand(x1 = 0, x2 = 1, size = 1, integer = False):
    #Get the random variable from 'rand.txt'
    rand = get_rand()
    if integer:
        x = np.zeros(size, dtype = np.int64)
    else:
        x = np.zeros(size)
    for i in range(size):
        #XOR-shift
        rand = rand ^ (rand >> 21)
        rand = rand ^ (rand << 37)
        rand = rand ^ (rand >> 4)
        #MCG
        rand = ((86729 * rand) % 2**32)
        if integer:
            x[i] = x1 + int(rand % (x2 - x1))
        else:
            x[i] = x1 + (x2 - x1) * rand / 2**32
    #Save the random variable to 'rand.txt'
    set_rand(rand)
    if size == 1:
        return x[0]
    else:
        return x

#Pearson correlation coefficient
def pearson(x, y):
    return sum((x - 0.5) * (y - 0.5)) / (np.sqrt(sum((x - 0.5)**2)) * np.sqrt((sum((y - 0.5)**2))))

#Hernquist profile
def hernquist(r):
    return M_dm / (2 * np.pi) * a / (r * (r + a)**3)

#Hernquist profile derivative
def hernquist_der(r):
    return a * M_dm / (2 * np.pi) * -(a + 4*r) / (r**2 * (r + a)**4)

#Hernquist 3D profile
def hernquist_3D(r):
    return hernquist(r) * 4 * np.pi * r**2

#Cumulative distribution function of the Hernquist profile
def cdf_hernquist(r):
    return -a * (2*r + a) / ((r + a)**2) + 1

#Inverted CDF
def cdf_hernquist_inv(y):
    return a * np.sqrt(y) / (1 - np.sqrt(y))

#Calculates the fraction of numbers in a distribution that are less than r
def fraction(dist, r):
    count = 0
    for i in dist:
        if i < r:
            count += 1
    return count / len(dist)

#Calculates the derivative of f with central differences
def central_diff(f, x, h):
    return(f(x + h) - f(x - h)) / (2 * h)

#Calculates the derivative of f with ridder's method
def ridders(f, x, h, d, m):
    diffs = np.zeros(m)
    #Define the initial approximations
    for i in range(m):
        diffs[i] = central_diff(f, x, h)
        h = h / d
    #Combine pairs of approximations to get a better final approximation
    for i in range(1, m):
        for j in np.flip(np.arange(m - i)):
            diffs[i + j] = (4**i * diffs[i + j] - diffs[i + j - 1]) / (4**i - 1)
    return diffs[m - 1]

#Calculates the root of f with the newton raphson method
def newton(f, f_der, x, max_err, it):
    i = 0
    c = x
    while i <= it:
        #Find points closer to the root by calculating the root of a first order Taylor series in the current point 
        c = c - f(c)/f_der(c)
        if np.abs(f(c)) < max_err:
            break
        i += 1
    return c, i

#2D Hernquist potential
def hernquist_pot(vec):
    x, y = vec[0], vec[1]
    return -G * M_dm / (np.sqrt((x - 1.3)**2 + 2 * (y - 4.2)**2) + a)

#Sorts an array x with selection sort (with option to sort f(x))
def selection_sort_f(x, f = lambda a : a):
    x = np.array(x)
    N = len(x)
    for i in range(0, N-1):
        #set the initial minimum
        imin = i
        #If the value of an element is equal to the previous element, there is no need to check the rest of the array since it already has the lowest possible value. 
        if i > 0 and f(x[i - 1]) == f(x[i]):
            continue
        #Check the rest of the array for values lower than the initial minimum and switch it with the lowest value. 
        for j in range(i+1, N):
            if f(x[j]) < f(x[imin]):
                imin = j
        if imin != i:
            x[[i,imin]] = x[[imin,i]]
    return x

#Moves toward the minimum of f with the downhill simplex method
#Terminates when a certain error values is reached or after a number of iterations
def downhill_simplex(f, x0, max_err, its):
    x0 = np.array(x0)
    dim = len(x0)
    
    #Create the first simplex
    x = np.zeros((dim + 1, dim))
    x[0] = x0
    for i in range(1,dim + 1):
        init = np.zeros(2)
        init[i-1] = 5
        x[i] = x0 + init
        
    #Determine the best point and calculate the centroid
    x = selection_sort_f(x, f)
    centroid = sum(x[:-1])/(dim)
    
    #Alter the simplex to find points closer to the minimum
    it = 0
    while 2*np.abs(f(x[-1]) - f(x[0]))/ np.abs(f(x[-1]) + f(x[0])) > max_err:
        #try reflection
        x_try = 2 * centroid - x[-1]
        if f(x_try) >= f(x[0]) and f(x_try) < f(x[-1]):
            x[-1] = x_try
        elif f(x_try) < f(x[0]):
            #try expansion
            x_exp = 2 * x_try - centroid
            if f(x_exp) < f(x_try):
                x[-1] = x_exp
            else:
                x[-1] = x_try
        else:
            #try contraction
            x_try = (centroid + x[-1])
            if f(x_try) < f(x[-1]):
                x[-1] = x_try
            else:
                #nothing else works, shrink
                for i in range(1, dim + 1):
                    x[i] = (centroid + x[i])/2
        #Determine the new best point and calculate the centroid
        x = selection_sort_f(x, f)
        centroid = sum(x[:-1])/(dim)
        it += 1
        if it > its:
            break
    return x[0]

#Creates a distribution from function f with rejection sampling
def reject_sample(f, a, b, size):
    x = []
    while len(x) < size:
        #Create random point in a rectangle that contains the function 
        Ux = my_rand(a, b)
        Uy = my_rand(0, 1)
        #If the point is below the function, add to the distribution
        if Uy <= f(Ux):
            x.append(Ux)
    return x
    
