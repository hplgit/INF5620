from numpy import *

def leapfrog(I, a, T, dt, first_step='ForwardEuler'):
    N = int(round(T/float(dt)))  # no of intervals
    u = zeros(N+1)
    t = linspace(0, T, N+1)

    u[0] = I
    if first_step == 'ForwardEuler':
        u[1] = u[0] - dt*a*u[0]
    elif first_step == 'CrankNicolson':
        u[1] = (1 - 0.5*a*dt)/(1 + 0.5*dt*a)*u[0]
    for n in range(1, N):
        u[n+1] = u[n-1] - 2*dt*a*u[n]
    return u, t

def leapfrog_filtered(I, a, T, dt, gamma=0.06):
    """
    Filtered version of Leapfrog. gamma=0.06 correspond to
    NCAR Climate Model.
    """
    N = int(round(T/float(dt)))  # no of intervals
    u = zeros(N+1)
    t = linspace(0, T, N+1)

    u[0] = I
    u[1] = (1 - 0.5*a*dt)/(1 + 0.5*dt*a)*u[0]  # Crank-Nicolson 1. step
    for n in range(1, N):
        u[n+1] = u[n-1] - 2*dt*a*u[n]
        u[n] = u[n] + gamma*(u[n-1] - 2*u[n] + u[n+1])
    return u, t

def Backward2ndOrder(I, a, T, dt):
    """2nd-order backward scheme."""
    N = int(round(T/float(dt)))  # no of intervals
    u = zeros(N+1)
    t = linspace(0, T, N+1)

    u[0] = I
    u[1] = (1 - 0.5*a*dt)/(1 + 0.5*dt*a)*u[0]  # Crank-Nicolson 1. step
    for n in range(1, N):
        u[n+1] = (4*u[n] - u[n-1])/(3 + 2*dt*a)
    return u, t

def Adams_Bashforth_3(I, a, T, dt):
    """
    Third-order Adams-Bashforth scheme (4 time levels), with
    Crank-Nicolson 2nd order scheme for the first two steps.
    """
    N = int(round(T/float(dt)))  # no of intervals
    u = zeros(N+1)
    t = linspace(0, T, N+1)

    u[0] = I
    u[1] = (1 - 0.5*a*dt)/(1 + 0.5*dt*a)*u[0]  # Crank-Nicolson 1. step
    u[2] = (1 - 0.5*a*dt)/(1 + 0.5*dt*a)*u[1]  # Crank-Nicolson 2. step
    for n in range(2, N):
        u[n+1] = u[n] - dt/12.*a*(23*u[n] - 16*u[n-1] + 5*u[n-2])
    return u, t
