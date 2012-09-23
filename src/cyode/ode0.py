"""
ODE integration restricted to scalar ODEs.
No use of arrays.
"""
def solver(f, I, dt, T, method):
    """
    Solve scalar ODE: 
    u'(t) = f(u,t), u(0)=I, 0 < t <= T
    method: numerical method to advance u one time step.
    dt: time step length.
    """
    N = int(round(float(T)/dt))
    u = I
    t = 0
    for n in xrange(N):  # may get memory error for large N
        u = method(u, t, f, dt)
        t += dt
    return u, t

def RK2(u, t, f, dt):
    """
    2nd-orderRunge-Kutta method for advancing u at time t
    to a new value at time t+dt. f is the right-hand side
    function f(u,t) in the equation u'=f.
    """
    K1 = dt*f(u, t)
    K2 = dt*f(u + 0.5*K1, t + 0.5*dt)
    unew = u + K2
    return unew
    
def problem1(u, t):
    """Right-hand side function f(u,t) for the ODE u'=-u+1."""
    return -u + 1

from math import exp

def problem2(u, t):
    return - u + exp(-2*t)





    
