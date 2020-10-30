#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from ex1_module import selection_sort_f, my_rand

#Constants for the number density profile
a_n = 2.4
b_n = 0.25
c_n = 1.6 
N_sat = 100
A = 256/(5*np.pi**(1.5))

#Finds a bracket near a and b in which a minimum of function f lies
def bracket_finder(f, a, b):
    #Make sure f(b) is lower than f(a)
    if f(b) > f(a):
        a, b = b, a
    i = 0
    #Define w as the golden ratio
    w = (1 + np.sqrt(5))/2
    #Define a new point c
    c = b + (b - a) * w
    #Try to find a bracket where the middle point is lower than the other two points
    while f(c) < f(b):
        #Calculate point d as the abscissa of the minimum of a parabola going through a, b and c
        d = ( -(c**2 * (f(a) - f(b)) + b**2 * (f(c) - f(a)) + a**2 * (f(b) - f(c))) /
         (2 * (c * (f(b) - f(a)) + b * (f(a) - f(c)) + a * (f(c) - f(b)))))
        #d is between b and c 
        if d > b and d < c:
            #bracket [b,d,c]
            if f(d) < f(c):
                a = b
                b = d
            #bracket [a,b,d]
            elif f(d) > f(b):
                c = d
            #no bracket found, try looking farther than c
            else:
                d = c + (c - b) * w
        #d must be beyond c
        else:
            #take another step if d is too far
            if np.abs(d - b) <= 100 * np.abs(c - b):
                d = c + (c - b) * w
        #redefine a, b and c by moving down the slope
        if f(c) < f(b):
            a = b
            b = c
            c = d       
        i += 1
    if b > a:
        return [a, b, c]
    else:
        return [c, b, a]

#Finds the minimum of function f in a bracket with the golden search method
#The golden search method decreases the width of the bracket (with the minimum still inside it)
#until its width is smaller than max_err.
def golden_search(f, bracket, max_err): 
    i = 0
    #initialize the bracket
    a,b,c = bracket[0], bracket[1], bracket[2]
    w = 2 - (1 + np.sqrt(5))/2
    while np.abs(c - a) > max_err:
        #Calculate a point d inside the biggest interval of [a,b] and [b,c]
        if np.abs(c - b) > np.abs(b - a):
            d = b + (c - b) * w
        else:
            d = a + (b - a) * w
        #Tighten towards d if f(d) < f(b)
        if f(d) < f(b):
            if d > a and d < b:
                c = b
                b = d
            else:
                a = b
                b = d
        #Otherwise tighten towards b
        else:
            if d > a and d < b:
                a = d
            else:
                c = d
        i += 1
        if i > 1000:
            raise Exception('infinite loop')
    #Return the interval which contains the minimum
    if f(d) < f(b):
        return d, i
    else:
        return b, i

#Number density satellite profile
def n(x):
    return A * N_sat * (x/b_n)**(a_n - 3) * np.exp(-(x/b_n)**c_n)

#N(x)dx gives the number of satellites in dx
def N(x):
    return n(x) * 4 * np.pi * x**2

#Satellite probability distribution (normalized with N_sat = 100)
def p(x):
    return N(x)/ N_sat

#The negative of N(x) (in order to create a function with a minimum instead of a maximum)
def N_neg(x):
    return -1 * N(x)

#Selects 'size' random elements from an array arr and returns them as a list
def rand_pick(arr, size):
    N = len(arr)
    if size > N:
        raise Exception('Cannot pick {} random elements from array of size {}'.format(size, N))
    #Copy the array
    dummy = list(arr)
    x = []
    for i in np.arange(size):
        #Make random index to pick a random value from the array, account for the fact that the dummy array gets smaller
        ind = my_rand(0, N-1-i, integer = True)
        #Remove used values from the dummy list, so they can't be used again.
        x.append(dummy.pop(ind))
    return x

#Gives the number of galaxies in with a radius less than 'r'
def within_r(arr, r):
    return len([gal for gal in arr if gal <= r])

#Returns a distribution of points of size 'size' based on a function f between a and b
#The function can be normalized with 'norm'
def reject_sample(f, a, b, size, norm = 1):
    x = []
    while len(x) < size: 
        Ux = my_rand(a, b)
        Uy = my_rand(0, 1)
        #Option to normalize the function
        if Uy <= f(Ux)*norm:
            x.append(Ux)
    return np.array(x)

#Returns the 'counts' in 'bins' from a 'distribution'
def histogram(distribution, bins):
    counts = np.zeros(len(bins) - 1, dtype = np.int64)
    for gal in distribution:
        bin_index = 0
        #While the galaxy exeeds the bin, go to the next bin
        while gal >= bins[bin_index]:
            bin_index += 1
        counts[bin_index - 1] += 1
    return counts

#Randomly divides the elements in array arr over sub_arrays of a certain size
def divide_in_arrays(arr, size):
    N = len(arr)
    if size > N:
        raise Exception('Cannot pick {} random elements from array of size {}'.format(size, N))
    #Make a copy of the array
    dummy = list(arr)
    number_of_arrays = int(N/size)
    x = np.zeros((number_of_arrays, size))
    for i in np.arange(number_of_arrays):
        for j in np.arange(size):
            #Make random index to pick a random value from the array, account for the fact that the dummy array gets smaller
            ind = my_rand(0, N-(size*i)-j, integer = True)
            #Remove used values from the dummy list, so they can't be used again.
            x[i,j] = dummy.pop(ind)
    return x

#Formula for Poisson distribution (taken from assignment 1)
def P(lam, k):
    #To avoid the overflow from dealing with large numbers calculate the Poisson distribution in log space
    #and then convert back to normal space
    #The factorial in log space is simply the summation of a list with the logarithm of all (positive) integers up to k
    return 10**(k*np.log10(lam) + np.log10(np.e) * (-lam) - sum(np.log10(np.arange(1, k+1))))

