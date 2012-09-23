"""Solve system of ODEs."""

from scitools.std import plot, figure

import numpy as np 
cimport numpy as np
cimport cython
ctypedef np.float64_t DT

cdef class Problem:
    # dudt is not available in Python space, only in C space
    #buffer types only allowed as local function variables
    #cdef np.ndarray[DT, ndim=1, negative_indices=False,
    #                mode='c'] dudt
    cdef np.ndarray dudt
    cdef np.ndarray u
    # remedy: send dudt to rhs as in F77/C
    
    def __init__(self):
        self.dudt = np.zeros(2)

    def rhs(self,
            np.ndarray[DT, ndim=1, negative_indices=False,
                       mode='c'] u,
            double t):
        return 0

cdef class Problem1(Problem):
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
RK2 = Method_RK2()
problem1 = Problem1()

@cython.boundscheck(False) # turn off bounds checking for this func.
def solver(Problem f, U0_, t_, ODEMethod method):
    # U0_ and t_ can be flexible objects
    cdef np.ndarray[DT, ndim=1, negative_indices=False,
                    mode='c'] t = np.asarray(t_)
    N = len(t_)-1
    if isinstance(U0_, (float,int)):
        U0_ = [U0_]  # wrap in list, which then will be array
    cdef np.ndarray[DT, ndim=1, negative_indices=False,
                    mode='c'] U0 = np.asarray(U0_)
    if not isinstance(f.rhs(U0,0), np.ndarray):
        raise TypeError('f (%s) must return numpy array' % f.__name__)

    cdef np.ndarray[DT, ndim=2, negative_indices=False,
                    mode='c'] u = np.zeros((N+1, len(U0)))
    u[0,:] = U0[:]
                    
    for n in range(N):
        u[n+1,:] = method.advance(u, n, t, f)
    return u, t




    
