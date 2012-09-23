"""Solve system of ODEs."""

import numpy as np
from scitools.std import plot, figure

def solver(f, I, t, method):
    t = np.asarray(t)
    N = len(t)-1
    if isinstance(I, (float,int)):
        I = [I]  # wrap in list, which then will be array
    I = np.asarray(I)
    if not isinstance(f(I,0), np.ndarray):
        raise TypeError('f (%s) must return numpy array' % f.__name__)
    u = np.zeros((N+1, len(I)))
    u[0] = I[:]

    for n in range(N):
        u[n+1] = method(u, n, t, f)
    return u, t

def RK2(u, n, t, f):
    dt = t[n+1] - t[n]
    K1 = dt*f(u[n], t[n])
    K2 = dt*f(u[n] + 0.5*K1, t[n] + 0.5*dt)
    unew = u[n] + K2
    return unew

def problem2(u, t):
    r = [u[1], -u[0]]
    return np.asarray(r)

# Improvement: send dudt array to f as optional argument for effciency
dudt = np.zeros(2)  # storage used by problem2

def problem3(u, t):
    dudt[0] = u[1]
    dudt[1] =-u[0]
    return dudt

class Problem1:
    def __init__(self):
        # Allocate an array for dudt for efficiency
        self.dudt = np.zeros(2)

    def __call__(self, u, t):
        self.dudt[0] = u[1]
        self.dudt[1] = -u[0]
        return self.dudt

def case(nperiods=4, showplot=False, ftype='class'):
    I = [1, 0]
    if ftype == 'class':
        f = Problem1()
    elif ftype == 'func2':
        f = problem2
    elif ftype == 'func3':
        f = problem3
    t0 = time.clock()
    time_points = np.linspace(0, nperiods*2*np.pi, nperiods*30+1)
    u, t = solver(f, I, time_points, RK2)
    t1 = time.clock()
    if showplot:
        plot(t[-200:], u[-200:,0])
        figure()
        plot(u[-200:,0], u[-200:,1])
    return t1-t0, time_points.size

if __name__ == '__main__':
    import time
    nperiods = 1000
    cpu, n = case(nperiods, showplot=True, ftype='class')
    print 'CPU-time, class version    :', cpu, 's', n
    cpu, n = case(nperiods, showplot=True, ftype='func1')
    print 'CPU-time, function version1:', cpu, 's', n
    cpu, n = case(nperiods, showplot=True, ftype='func2')
    print 'CPU-time, function version2:', cpu, 's', n






