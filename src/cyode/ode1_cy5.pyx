"""
Cython code restricted to scalar ODEs.
Variables are declared with types.
Functions as arguments are represented by classes and instances.
Numpy arrays are declared with 1) fixed number of dimensions, 
2) element type, 3) negative indices turned off, 4) bounds checking
off, and 5) contiguous memory.
Compared to ode1_cy4, here double is replaced by np.float_t. (No
significant change in efficiency.)
"""
import numpy as np 
cimport numpy as np
cimport cython

cdef class Problem:
    cpdef np.float_t rhs(self, np.float_t u, np.float_t t):
        return 0

cdef class Problem1(Problem):
    cpdef np.float_t rhs(self, np.float_t u, np.float_t t):
        return -u +1  # u = 1-exp(-t)

cdef extern from "math.h":
    np.float_t exp(np.float_t)

cdef class Problem2(Problem):
    cpdef np.float_t rhs(self, np.float_t u, np.float_t t):
        return - u + exp(-2*t)

ctypedef np.float64_t DT

cdef class ODEMethod:
    cpdef advance(self, np.float_t u_1, int n, np.float_t t_1,
                  np.float_t dt, Problem p):
        return 0

cdef class Method_RK2(ODEMethod):
    cpdef advance(self, np.float_t u_1, int n, np.float_t t_1,
                  np.float_t dt, Problem p):
        cdef np.float_t K1, K2, unew
        K1 = dt*p.rhs(u_1, t_1)
        K2 = dt*p.rhs(u_1 + 0.5*K1, t_1 + 0.5*dt)
        unew = u_1 + K2
        return unew
    
# Create names compatible with ode1.py
RK2 = Method_RK2()
problem1 = Problem1()
problem2 = Problem2()

@cython.boundscheck(False) # turn off bounds checking for this func.
def solver(Problem f, 
           np.float_t I, 
           np.ndarray[DT, ndim=1, negative_indices=False, mode='c'] t, 
           ODEMethod method):
    cdef int N = len(t)-1
    #cdef np.ndarray[DT, ndim=1, negative_indices=False, mode='c'] u = np.zeros(N+1, dtype=np.float_t)
    #Cython does not like type specification via dtype when the buffer
    #declares the type
    cdef np.ndarray[DT, ndim=1, negative_indices=False, 
                    mode='c'] u = np.zeros(N+1)
    u[0] = I   
             
    cdef int n
    for n in range(N):
        u[n+1] = method.advance(u[n], n, t[n], t[n+1]-t[n], f)
    return u, t

