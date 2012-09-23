"""
Cython code restricted to scalar ODEs.
Variables are declared with types.
Functions as arguments are represented by classes and instances.
"""
# Note: need both numpy imports!
import numpy as np
cimport numpy as np

cdef class Problem:
    cpdef double rhs(self, double u, double t):
        return 0

cdef class Problem1(Problem):
    cpdef double rhs(self, double u, double t):
        return -u +1  # u = 1-exp(-t)

from math import exp

cdef class Problem2(Problem):
    cpdef double rhs(self, double u, double t):
        return - u + exp(-2*t)

cdef class ODEMethod:
    cpdef double advance(self, np.ndarray u, int n,
                         np.ndarray t, Problem p):
        return 0

cdef class Method_RK2(ODEMethod):
    cpdef double advance(self, np.ndarray u, int n,
                         np.ndarray t, Problem p):
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

cpdef solver(Problem f, double U0, np.ndarray t, ODEMethod method):
    cdef int N = len(t)-1
    cdef np.ndarray u = np.zeros(N+1, dtype=np.float)
    u[0] = U0

    cdef int n
    for n in range(N):
        u[n+1] = method.advance(u, n, t, f)
    return u, t

