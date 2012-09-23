"""
ODE integration restricted to scalar ODEs.
No use of arrays.
Just plain Python code compiled as Cython.
"""

def solver(f, U0, dt, t_end, method):
    N = int(round(float(t_end)/dt))
    u = U0  # previous time step
    t = 0
    for n in xrange(N):
        u = method(u, t, f, dt)
        t += dt
    return u, t

def RK2(u, t, f, dt):
    K1 = dt*f(u, t)
    K2 = dt*f(u + 0.5*K1, t + 0.5*dt)
    unew = u + K2
    return unew
    
def problem1(u, t):
    return -u +1  # u = 1-exp(-t)

from math import exp

def problem2(u, t):
    return - u + exp(-2*t)





    
