"""Unified scalar/system ODE implementation."""

import numpy as np
from scitools.std import plot, figure

def solver(f_, U0, t, method):
    f = lambda u, t: np.asarray(f_(u, t))
    t = np.asarray(t)
    N = len(t)-1
    if isinstance(U0, (list,tuple,np.ndarray)):
        # System of ODEs
        u = np.zeros((N+1, len(U0)))
        u[0] = np.asarray(U0)
    else:
        # Scalar ODE
        u = np.zeros(N+1)
        u[0] = U0
                    
    for n in range(N):
        u[n+1] = method(u, n, t, f)
    return u, t

def RK2(u, n, t, f):
    """2nd-orderRunge-Kutta method."""
    dt = t[n+1] - t[n]
    K1 = dt*f(u[n], t[n])
    K2 = dt*f(u[n] + 0.5*K1, t[n] + 0.5*dt)
    unew = u[n] + K2
    return unew
    

def problem(u, t):
    return [u[1], -u[0]]

class ProblemOpt:
    """RHS of ODE system. Allocates an array for dudt for efficiency."""
    def __init__(self):
        self.dudt = np.zeros(2)
        
    def __call__(self, u, t):
        self.dudt[0] = u[1]
        self.dudt[1] = -u[0]
        return self.dudt
    
def case(nperiods=4, showplot=False, ftype='class'):
    U0 = [1, 0]
    if ftype == 'class':
        f = ProblemOpt()
    else:
        f = problem
    t0 = time.clock()
    u, t = RK2(f, U0, t=np.linspace(0, nperiods*np.pi, nperiods*30+1))
    t1 = time.clock()
    if showplot:
        plot(t[-200:], u[-200:,0])
        figure()
        plot(u[-200:,0], u[-200:,1])
    return t1-t0

if __name__ == '__main__':
    import time
    cpu = case(1000, showplot=True, ftype='class')
    print 'CPU-time, class version   :', cpu
    cpu = case(1000, showplot=True, ftype='func')
    print 'CPU-time, function version:', cpu





    
