"""ODE integration restricted to scalar ODEs."""
from scitools.std import plot, figure
import numpy as np

def solver(f, I, t, method):
    t = np.asarray(t)
    N = len(t)-1
    u = np.zeros(N+1)
    u[0] = I

    for n in range(N):
        u[n+1] = method(u, n, t, f)
    return u, t

def RK2(u, n, t, f):
    dt = t[n+1] - t[n]
    K1 = dt*f(u[n], t[n])
    K2 = dt*f(u[n] + 0.5*K1, t[n] + 0.5*dt)
    unew = u[n] + K2
    return unew

def problem1(u, t):
    return -u + 1

from math import exp

def problem2(u, t):
    return -u + exp(-2*t)






