#!/usr/bin/env python3
import numpy as np

#Importing the data
w_gs = np.loadtxt("wgs.dat",dtype=np.float32)
w_ss = np.loadtxt("wss.dat",dtype=np.float32)

#Solves a matrix equation Ax = b with LU-decomostion for x, given A and b.
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

