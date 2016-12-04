'es# -*- coding: utf-8 -*-
"""
Created on Sun Nov 13 20:19:05 2016
#from http://hplgit.github.io/teamods/cyode/cyode-sphinx/main_cyode.html
#http://hplgit.github.io/num-methods-for-PDEs/doc/pub/index.html
@author: Hans Petter Langtangen ()
"""
import os
import numpy as np
cimport numpy as np
cimport cython

ctypedef np.float64_t DT




#
#

cdef class Problem:             #cdef class will define a "class of object"
    cdef np.ndarray dudt        #the class will contain a ndarray called dudt

    def __init__(self):         #these 2 lines will assign a value to dudt  (self refers to an object of the "Problem" class. Let say later on in the code, you create have a line ")
        self.dudt = np.zeros(2)

    def rhs(self,
            np.ndarray[DT, ndim=1, negative_indices=False,
                       mode='c'] u,
            double t):
        return 0

cdef class Problem1(Problem):  #we define a "subclass" of the class Problem
    def rhs(self,
            np.ndarray[DT, ndim=1, negative_indices=False,
                       mode='c'] u,
            double t):
        self.dudt[0] = u[1]
        self.dudt[1] = -u[0]
        return self.dudt


cdef class ODEMethod:
    def advance(self,
                np.ndarray[DT, ndim=2, negative_indices=False,
                           mode='c'] u,
                int n,
                np.ndarray[DT, ndim=1, negative_indices=False,
                           mode='c'] t,
                Problem p):
        return 0

@cython.boundscheck(False)
cdef class Method_RK2(ODEMethod):
    def advance(self,
                np.ndarray[DT, ndim=2, negative_indices=False,
                           mode='c'] u,
                int n,
                np.ndarray[DT, ndim=1, negative_indices=False,
                           mode='c'] t,
                Problem p):
        cdef np.ndarray[DT, ndim=1, negative_indices=False,
                        mode='c'] K1, K2, unew
        cdef double dt
        cdef np.ndarray[DT, ndim=1, negative_indices=False,
                        mode='c'] un = u[n,:]
        dt = t[n+1] - t[n]
        K1 = dt*p.rhs(un, t[n])
        K2 = dt*p.rhs(un + 0.5*K1, t[n] + 0.5*dt)
        unew = u[n,:] + K2
        return unew

# Create names compatible with ode2.py
RK2 = Method_RK2()               #this means that RK2 is an object of the class Method_RK2
problem1 = Problem1()            #this means that problem1 is an object of the class Problem1

@cython.boundscheck(False) # turn off bounds checking for this func.
def solver(Problem f, I_, t_, ODEMethod method):
    # I_ and t_ can be flexible objects
    cdef np.ndarray[DT, ndim=1, negative_indices=False,
                    mode='c'] t = np.asarray(t_)
    N = len(t_)-1
    if isinstance(I_, (float,int)):
        I_ = [I_]  # wrap in list, which then will be array
    cdef np.ndarray[DT, ndim=1, negative_indices=False,
                    mode='c'] I = np.asarray(I_)
    if not isinstance(f.rhs(I,0), np.ndarray):
        raise TypeError('f (%s) must return numpy array' %
                        f.__name__)

    cdef np.ndarray[DT, ndim=2, negative_indices=False,
                    mode='c'] u = np.zeros((N+1, len(I)))
    u[0,:] = I[:]

    for n in range(N):
        u[n+1,:] = method.advance(u, n, t, f)
    return u, t


    


