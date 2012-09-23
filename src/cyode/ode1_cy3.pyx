"""
Cython code restricted to scalar ODEs.
Variables are declared with types.
Functions as arguments are represented by classes and instances.
Numpy arrays are declared with fixed dimensions and type.
The RK2 function avoids use of arrays and have scalar arguments
only (thus allowing cpdef and faster calls).
"""
import numpy as np   # note: need both imports!
cimport numpy as np

cdef class Problem:
    cpdef double rhs(self, double u, double t):
        return 0

cdef class Problem1(Problem):
    cpdef double rhs(self, double u, double t):
        return -u +1  # u = 1-exp(-t)

#from math import exp
cdef extern from "math.h":
    double exp(double)

#or
#from libc.math cimport exp  # may need explicit -lm linking
#see http://docs.cython.org/src/tutorial/external.html

cdef class Problem2(Problem):
    cpdef double rhs(self, double u, double t):
        return - u + exp(-2*t)

# NOTE: need def, not cpdef, for functions with array arguments
# and [] buffer notation.
# This means that def functions with arrays are not called very
# efficiently, and the RK2.advance function, which basically
# works with a single array element should be implemented alternatively
# via doubles only.

cdef class ODEMethod:
    cpdef advance(self, double u_1, int n, double t_1,
                  double dt, Problem p):
        return 0

cdef class Method_RK2(ODEMethod):
    cpdef advance(self, double u_1, int n, double t_1,
                  double dt, Problem p):
        cdef double K1, K2, unew
        K1 = dt*p.rhs(u_1, t_1)
        K2 = dt*p.rhs(u_1 + 0.5*K1, t_1 + 0.5*dt)
        unew = u_1 + K2
        return unew
    
# Create names compatible with ode0.py
RK2 = Method_RK2()
problem1 = Problem1()
problem2 = Problem2()


def solver(Problem f, double I, 
           np.ndarray[np.float_t, ndim=1] t, 
           ODEMethod method):
    cdef int N = len(t)-1
    #cdef np.ndarray[np.float_t, ndim=1] u = np.zeros(N+1, dtype=np.float_t)
    #Cython does not like type specification via dtype when the buffer
    #declares the type
    cdef np.ndarray[np.float_t, ndim=1] u = np.zeros(N+1)
    u[0] = I   
             
    cdef int n
    for n in range(N):
        u[n+1] = method.advance(u[n], n, t[n], t[n+1]-t[n], f)
    return u, t

