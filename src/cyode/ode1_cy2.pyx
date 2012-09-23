"""
Cython code restricted to scalar ODEs.
Variables are declared with types.
Functions as arguments are represented by classes and instances.
Numpy arrays are declared with fixed dimensions and type.
"""
import numpy as np   # note: need both imports!
cimport numpy as np

cdef class Problem:
    cpdef double rhs(self, double u, double t):
        return 0

cdef class Problem1(Problem):
    cpdef double rhs(self, double u, double t):
        return -u +1  # u = 1-exp(-t)

cdef extern from "math.h":
    double exp(double)

cdef class Problem2(Problem):
    cpdef double rhs(self, double u, double t):
        return - u + exp(-2*t)

# NOTE: need def, not cpdef, for functions with array arguments
# and [] buffer notation.
# Common error message: "Expected ']'"

cdef class ODEMethod:
    def advance(self, 
                np.ndarray[np.float_t, ndim=1] u,
                int n, 
                np.ndarray[np.float_t, ndim=1] t, 
                Problem p):
        return 0

cdef class Method_RK2(ODEMethod):
    def advance(self, 
                np.ndarray[np.float_t, ndim=1] u, 
                int n,
                np.ndarray[np.float_t, ndim=1] t, 
                Problem p):
        """2nd-orderRunge-Kutta method."""
        cdef double K1, K2, unew, dt
        dt = t[n+1] - t[n]
        K1 = dt*p.rhs(u[n], t[n])
        K2 = dt*p.rhs(u[n] + 0.5*K1, t[n] + 0.5*dt)
        unew = u[n] + K2
        return unew
    
# Create names compatible with ode0.py
RK2 = Method_RK2()
problem1 = Problem1()
problem2 = Problem2()


def solver(Problem f, double U0, 
           np.ndarray[np.float_t, ndim=1] t, 
           ODEMethod method):
    cdef int N = len(t)-1
    #cdef np.ndarray[np.float_t, ndim=1] u = np.zeros(N+1, dtype=np.float_t)
    #Cython does not like type specification via dtype when the buffer
    #declares the type
    cdef np.ndarray[np.float_t, ndim=1] u = np.zeros(N+1)
    u[0] = U0   
             
    cdef int n
    for n in range(N):
        u[n+1] = method.advance(u, n, t, f)
    return u, t

