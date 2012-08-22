from numpy import *
from matplotlib.pyplot import *

def theta_rule(I, a, b, T, dt, theta):
    """
    Solve u'=-a(t)*u + b(t), u(0)=I,
    for t in (0,T] with steps of dt.
    a and b are Python functions of t.
    """
    N = int(round(T/float(dt)))  # no of intervals
    u = zeros(N+1)
    t = linspace(0, T, N+1)

    u[0] = I
    for n in range(0, N):
        u[n+1] = ((1 - dt*(1-theta)*a(t[n]))*u[n] + \
                  dt*(theta*b(t[n+1]) + (1-theta)*b(t[n])))/\
                  (1 + dt*theta*a(t[n+1]))
    return u, t

def verify_linear_solution():
    def exact_solution(t):
        return c*t + I

    def a(t):
        return t**0.5  # can be arbitrary

    def b(t):
        return c + a(t)*exact_solution(t)

    theta = 0; I = 0.1; dt = 0.1; c = -0.5
    T = 4
    N = int(T/dt)  # no of steps
    u, t = theta_rule(I=I, a=a, b=b, T=N*dt, dt=dt, theta=theta)
    u_e = array([exact_solution(tn) for tn in t])
    difference = abs(u_e - u).max()  # max deviation
    tol = 1E-15  # tolerance for comparing floats
    success = difference <= tol
    return success

if __name__ == '__main__':
    print verify_linear_solution()


