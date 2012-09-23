"""
ODE integration restricted to scalar ODEs.
No use of arrays.
Cython version with declaration of variables.
"""

cpdef solver(f, double U0, double dt, double t_end, method) except *:
    cdef int N = int(round(float(t_end)/dt))
    cdef double u = U0  # previous time step
    cdef double t = 0
    cdef int n
    for n in xrange(N):
        u = method(u, t, f, dt)
        t += dt
    return u, t

cpdef double RK2(double u, double t, f, double dt) except *:
    cdef double K1, K2, unew
    K1 = dt*f(u, t)
    K2 = dt*f(u + 0.5*K1, t + 0.5*dt)
    unew = u + K2
    return unew
    
cpdef double problem1(double u, double t) except *:
    return -u +1  # u = 1-exp(-t)






    
