"""
ODE integration restricted to scalar ODEs.
No use of arrays.
Cython version with declaration of variables.
No use a functions transferred as arguments.
No use of except * in cpdef functions (factor 2.5).
"""

cpdef solver(f, double U0, double dt, double t_end, method):
    cdef int N = int(round(float(t_end)/dt))
    cdef double u = U0  # previous time step
    cdef double t = 0
    cdef int n
    for n in xrange(N):
        #u = method(u, t, f, dt)
        u = RK2(u, t, f, dt)
        t += dt
    return u, t

cpdef double RK2(double u, double t, f, double dt):
    cdef double K1, K2, unew
    #K1 = dt*f(u, t)
    #K2 = dt*f(u + 0.5*K1, t + 0.5*dt)
    K1 = dt*problem1(u, t)
    K2 = dt*problem1(u + 0.5*K1, t + 0.5*dt)
    unew = u + K2
    return unew
    
cpdef double problem1(double u, double t):
    return -u + 1

from math import exp

cpdef double problem2(double u, double t):
    return - u + exp(-2*t)





    
