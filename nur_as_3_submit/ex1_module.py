#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from scipy.special import gammainc

"""Functions"""

#N(x)dx gives the number of satellites in dx
def N(x, par, N_sat):
    a, b, c = par[0], par[1], par[2]
    return 4 * np.pi * N_sat * (x)**(a - 1) / (b)**(a - 3) * np.exp(-(x/b)**c)

#Derivative of N(x) with respect to a, b or c
def N_der(x, par, der_par, N_sat):
    a, b, c = par[0], par[1], par[2]
    
    #partial derivative with respect to a
    if der_par == 0: 
        return 4 * np.pi * N_sat * np.log(x/b) * (x)**(a - 1) / (b)**(a - 3) * np.exp(-(x/b)**c)
    
    #partial derivative with respect to b
    elif der_par == 1: 
        return 4 * np.pi * N_sat * (x)**(a - 1) * (3 - a) * (b)**(2 - a) * np.exp(-(x/b)**c) * c * x**c * b**(-c - 1)
    
    #partial derivative with respect to c
    elif der_par == 2: 
        return 4 * np.pi * N_sat * (x)**(a - 1) / (b)**(a - 3) * np.exp(-(x/b)**c) * np.log(b/x) * (x/b)**c
    else:
        raise ValueError('Function has only three adjustable parameters')

#Chi squared 
def chi_sq(data, model, var):
    return sum((data - model)**2 / var)

#First derivative of Chi squared
def chi_sq_1st_der(data, model, var, model_der):
    return -2 * sum((data - model) / var * model_der)

#Second derivative of Chi squared
def chi_sq_2nd_der(var, model_der_1, model_der_2):
    return 2 * sum(1/var * model_der_1 * model_der_2)

#Chi squared cumulative distribution function
def chi2_cdf(chi2, df):
    return gammainc(df/2, chi2/2)

#Poisson minus log likelihood
def min_ln_lik_poisson(data, model):
    return -1*sum(data*np.log(model) - model)

#Test statistic G
def G_test(data, model):
    #First remove empty bins
    model_ = np.delete(model, np.where(data == 0))
    data_ = np.delete(data, np.where(data == 0))
    #Degrees of freedom (the minus 4 is because of the three free parameters the last bin being constrained by the total number of satellites)
    df = len(model_) - 4
    return 2 * sum(data_ * np.log(data_/model_)), df

#Significance Q of a test statistic (with underlying Gaussian distribution)
def Q_value(test_statistic, df):
    return 1 - chi2_cdf(test_statistic, df)

#Returns the 'counts' in 'bins' from a 'distribution'
def histogram(distribution, bins):
    counts = np.zeros(len(bins) - 1, dtype = np.int64)
    for gal in distribution:
        bin_index = 0
        while gal >= bins[bin_index]:
            bin_index += 1
        counts[bin_index - 1] += 1
    return counts

#Retrieves the data from the files and creates useful arrays
def sats_list(filename):
    with open(filename, 'r') as f:
        haloes = [line for line in f.read().split('\n')]

    n_haloes = int(haloes[3])
    r = []
    theta = []
    phi = []
    haloes_with_sats = []

    has_sats = False
    for halo in haloes[4:]:
        if halo != '#' and halo != '':
            sat = np.array([[float(x) for x in halo.split('   ')]], dtype=np.float32)
            r.append(sat[0][0])
            if has_sats:
                sats = np.vstack((sats, sat))
            if not has_sats:
                sats = sat             
                has_sats = True
        if has_sats and (halo == '#' or halo == ''):
            haloes_with_sats.append(sats)
            has_sats = False
    return haloes_with_sats, np.array(r), n_haloes

#Romberg integration method
def romberg(f, a, b, m):
    h = b - a
    r = np.zeros(m)
    r[0] = h * 0.5 * (f(a) + f(b))
    for i in range(1, m):
        h = 0.5 * h
        r[i] = 0.5 * r[i - 1] + h * sum([f(a + (2*k - 1) * h) for k in range(1, 2**(i-1)+1)])
    for i in range(1, m):
        for j in np.flip(np.arange(m - i)):
            r[i + j] = (4**i * r[i + j] - r[i + j - 1]) / (4**i - 1)
    return r[m - 1]

#Calculates the value of the model for a set of bins
def model_vals(par, N_sat, n_haloes):
    
    #Define new function that only takes one variable for the romberg integration
    def new_N(x):
        return N(x, par, N_sat)
    
    #Boundaries for normalization
    x1 = 1e-4
    x2 = 5
    
    #Order of the romberg method, the normalization requires a high precision for the G-test  
    d1 = 13
    d2 = 5
    
    #Calculate the normalization
    A = N_sat/romberg(new_N, x1, x2, d1)
    
    model_size = len(bins) - 1
    model = np.zeros(model_size)
    for i in range(model_size):
        model[i] = A * romberg(new_N, bins[i], bins[i+1], d2)
    return model*n_haloes

#Calculates the derivative of the model for a set of bins
def model_der_vals(par, der_par, N_sat, n_haloes):
    
    #Define new functions that only take one variable for the romberg integration
    def new_N(x):
        return N(x, par, N_sat)
    def new_N_der(x):
        return N_der(x, par, der_par, N_sat)
    
    #Boundaries for normalization
    x1 = 1e-4
    x2 = 5
    
    #Order of the romberg method, 
    d1 = 13
    d2 = 5
    
    #Calculate the normalization
    A = N_sat/romberg(new_N, x1, x2, d1)
    
    #Initiate model
    model_size = len(bins) - 1
    model = np.zeros(model_size)
    
    #Calculate model
    for i in range(model_size):
        model[i] = A * romberg(new_N_der, bins[i], bins[i+1], d2)
    return model*n_haloes

def LU_solve(A, b):
    X, Y = A.shape
    
    #Create the upper and lower matrices
    L = np.zeros((X,Y),dtype=np.float32)
    U = np.zeros((X,Y),dtype=np.float32)
    
    for i in range(X):
        L[i,i] = 1

    for j in range(Y):
        U[0,j] = A[0,j]
        for i in range(j+1):
            U[i,j] = A[i,j] - sum([L[i,k]*U[k,j] for k in range(i)])
        for i in range(j+1, Y):
            L[i,j] = 1 / U[j,j] * (A[i,j] - sum([L[i,k] * U[k,j] for k in range(j)]))

    #Solve the matrix equation (LU)x = b
    y = np.zeros(Y,dtype=np.float32)
    x = np.zeros(X,dtype=np.float32)
    
    y[0] = b[0]/L[0,0]
    for i in range(1, X):
        y[i] = 1/L[i,i] * (b[i] - sum([L[i,j]*y[j] for j in np.arange(0, i)]))
    
    x[-1] = y[-1]/U[-1,-1]
    for i in np.flip(range(0, X)):
        x[i] = 1/U[i,i] * (y[i] - sum([U[i,j]*x[j] for j in np.arange(i+1, X)]))

    return L, U, x

def levenberg_marquardt(data, par0, min_check, max_iter, N_sat, n_haloes):
    """
    Function that minimizes the chi squared by adjusting the parameters of the model. 
    Parameters are updated via a Levenberg-Marquardt method.
    
    data: array with the y values of the data
    par0: array or list with the initial parameters in the form [a,b,c]
    min_chi2: minimum check, below which the algorithm terminates
    max_iter: maximum number of iterations 
    """
    
    par = np.array(par0)
    dim = len(par) #dimensions
    par_sets = []
    
    
    #Calculate the initial model values, variance and chi squared
    model = model_vals(par, N_sat, n_haloes)
    var = model
    chi2 = chi_sq(data, model, var)
    
    #Initialize variables for the levenberg marquardt method
    alpha = np.zeros((dim, dim))
    beta = np.zeros(dim)
    lam = 1e-3
    w = 10
    
    #Iteration variable
    i = 0
    no_new_par_counter = 0
    
    while i < max_iter:
        model_ders = []
        
        #Calculate the beta vector
        for k in range(dim):
            model_ders.append(model_der_vals(par, k, N_sat, n_haloes))
            beta[k] = -1/2 * chi_sq_1st_der(data, model, var, model_ders[k])
        
        #Calculate the alpha prime matrix
        for k in range(dim):
            for l in range(dim):
                alpha[k, l] = 1/2 * chi_sq_2nd_der(var, model_ders[k], model_ders[l])
                
                #
                if k == l:
                    alpha[k, l] = (1 + lam) * alpha[k, l]
        
        #Solve the matrix equation to calculate the change in the parameters
        dpar = LU_solve(alpha, beta)[2]
        
        #New parameters and associated variables
        new_par = par + dpar
        new_model = model_vals(new_par, N_sat, n_haloes)
        new_var = new_model
        new_chi2 = chi_sq(data, new_model, new_var)
        
        #If the new chi squared is greater, then move towards steepest descent
        if new_chi2 >= chi2:
            lam = lam * w
            no_new_par_counter += 1
        
        #If the new chi squared is smaller, then move towards Newton's method 
        #The parameters have improved so we update the parameters and associated variables
        else:
            no_new_par_counter = 0
            lam = lam / w
            par = new_par
            par_sets.append(par)
            model = new_model
            var = new_var
            chi2 = new_chi2
        
        #Breakchecks: the last 4 new(!) parameter sets are too close together 
        #or no new parameter set was found for 20 iterations
        i += 1
        if len(par_sets) >= 5:
            check = sum([sum(np.abs(np.array(par_sets)[-i] - np.array(par_sets)[-i-1])) for i in range(1, 5)])
            if check < min_check:
                break
        if no_new_par_counter > 20:
            break
    return par, chi2

#Moves toward the minimum of f with the downhill simplex method
#Terminates when a certain error values is reached or after a number of iterations
def downhill_simplex(f, x0, min_err, its):
    x0 = np.array(x0)
    dim = len(x0)
    best_points = []
    
    #Create the first simplex
    x = np.zeros((dim + 1, dim))
    x[0] = x0
    for i in range(1,dim + 1):
        init = np.zeros(dim)
        init[i-1] = 0.1*x0[i-1]
        x[i] = x0 + init
        
    #Determine the best point and calculate the centroid
    x = selection_sort_f(x, f)
    centroid = sum(x[:-1])/(dim)
    best_val = f(x[0])
    worst_val = f(x[-1])
    
    #Alter the simplex to find points closer to the minimum
    it = 0
    while 2*np.abs(worst_val - best_val)/ np.abs(worst_val + best_val) > min_err:
        old_best_point = x[0]
        
        #Make a new point by reflecting
        x_try = 2 * centroid - x[-1]
        
        #Evaluate the best, worst and new point
        best_val = f(x[0])
        worst_val = f(x[-1])
        try_val = f(x_try)
        
        #New point is better than the worst point, but not the best, accept it
        if try_val >= best_val and try_val < worst_val:
            x[-1] = x_try
        
        #New point is better than the best point, expand in same direciton
        elif try_val < best_val:
            x_exp = 2 * x_try - centroid
            
            #Expanded point is even better, accept expanded point
            if f(x_exp) < try_val:
                x[-1] = x_exp
            #Expanded point is worse, accept reflected point
            else:
                x[-1] = x_try
        
        #New point is worse than the worst point
        else:
            #Try contracting instead of reflecting
            x_try = centroid + x[-1]
            #Contraction point is better than worst point, accept it
            if f(x_try) < worst_val:
                x[-1] = x_try
            #No better points found, shrink towards the best point
            else:
                for i in range(1, dim + 1):
                    x[i] = (centroid + x[i])/2
                    
        #Reorder the points from best to worst and calculate new centroid
        x = selection_sort_f(x, f)
        centroid = sum(x[:-1])/(dim)
        
        #New evaluations
        best_val = f(x[0])
        worst_val = f(x[-1])
        
        #Breakchecks: max iterations reached or the best point is the same after 10 iterations
        it += 1
        best_points.append(x[0])
        if it > its:
            break
        if it > 11 and (best_points[-1][0] == best_points[-10][0]
                    and best_points[-1][1] == best_points[-10][1]
                    and best_points[-1][2] == best_points[-10][2]):
            break
    return x[0]

#Sorts an array x with selection sort (with option to sort f(x))
def selection_sort_f(x, f = lambda a : a):
    x = np.array(x)
    N = len(x)
    for i in range(0, N-1):
        imin = i
        if i > 0 and f(x[i - 1]) == f(x[i]):
            continue
        for j in range(i+1, N):
            if f(x[j]) <= f(x[imin]):
                imin = j
        if imin != i:
            x[[i,imin]] = x[[imin,i]]
    return x
    
"""Import data from files"""

#Create useful arrays from the files

filenames = sorted([file for file in glob.glob('satgals_m??.txt')])

#Total number of haloes (per file)
n_haloes_list = []

#The radius of all satellites (per file)
r_list = [] 

haloes_with_sats_list = []
                            
for filename in filenames:
    haloes_with_sats, r, n_haloes = sats_list(filename)
    n_haloes_list.append(n_haloes)
    haloes_with_sats_list.append(haloes_with_sats)
    r_list.append(r)

#Total number of satellites (per file)
tot_sats = [sum([len(halo) for halo in haloes_with_sats_list[i]]) for i in range(5)]

#Average number of satellites per halo (per file)
N_sats = [tot_sats[i]/n_haloes_list[i] for i in range(5)]

#Bins
n_bins = 20
bins = np.logspace(-4,np.log10(5), n_bins + 1)

#Calculate the counts for each bin for each file
counts = np.zeros((5, n_bins))
for m in range(5):
    counts[m] = np.array(histogram(np.array(r_list[m]), bins), dtype = np.float32)
