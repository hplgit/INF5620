import numpy as np

def FE_logistic(u0, dt, N):
    u = np.zeros(N+1)
    u[0] = u0
    for n in range(N):
        u[n+1] = u[n] + dt*(u[n] - u[n]**2)
    return u

def quadratic_roots(a, b, c):
    delta = b**2 - 4*a*c
    r2 = (-b + sqrt(delta))/float(2*a)
    r1 = (-b - sqrt(delta))/float(2*a)
    return r1, r2

def BE_logistic(u0, dt, Nt, choice='r2', eps_r=1E-3, omega=1):
    print 'BE_logistic', choice
    u = np.zeros(Nt+1)
    u[0] = u0
    for n in range(1, Nt+1):
        a = dt
        b = 1 - dt
        c = -u[n-1]
        if choice in ('r1', 'r2'):
            r1, r2 = quadratic_roots(a, b, c)
            u[n] = r1 if choice == 'r1' else r2
        elif choice == 'Picard':
            def F(u):
                return a*u**2 + b*u + c

            u_ = u[n-1]
            k = 0
            while abs(F(u_)) > eps_r:
                u_ = omega*(-c/(a*u_ + b)) + (1-omega)*u_
                k += 1
            print k, 'Picard iterations', F(u_)
            u[n] = u_
        elif choice == 'Newton':
            def F(u):
                return a*u**2 + b*u + c

            def dF(u):
                return 2*a*u + b

            u_ = u[n-1]
            k = 0
            while abs(F(u_)) > eps_r:
                u_ = u_ - F(u_)/dF(u_)
                k += 1
            print k, 'Newton iterations', F(u_)
            u[n] = u_
    return u

def CN_logistic(u0, dt, Nt):
    u = np.zeros(Nt+1)
    u[0] = u0
    for n in range(0, Nt):
        u[n+1] = (1 + 0.5*dt)/(1 + dt*u[n] - 0.5*dt)*u[n]
    return u

from scitools.std import *

def quadratic_root_goes_to_infinity():
    """
    Verify that one of the roots in the quadratic equation
    goes to infinity.
    """
    for dt in 1E-7, 1E-12, 1E-16:
        a = dt
        b = 1 - dt
        c = -0.1
        print dt, quadratic_roots(a, b, c)

import sympy as sp
dt, u_1, u = sp.symbols('dt u_1 u')
r1, r2 = sp.solve(dt*u**2 + (1-dt)*u - u_1, u)
print r1
print r2
print r1.series(dt, 0, 2)
print r2.series(dt, 0, 2)

dt = 0.8
N = 10
u_FE = FE_logistic(0.1, dt, N)
u_BE1 = BE_logistic(0.1, dt, N, 'r1')
u_BE2 = BE_logistic(0.1, dt, N, 'r2')
u_BE3 = BE_logistic(0.1, dt, N, 'Picard')
u_BE4 = BE_logistic(0.1, dt, N, 'Newton')
u_CN = CN_logistic(0.1, dt, N)

t = np.linspace(0, dt*N, N+1)
plot(t, u_FE, t, u_BE2, t, u_BE3, t, u_BE4, t, u_CN,
     legend=['FE', 'BE exact', 'BE Picard', 'BE newton', 'CN gm'])
savefig('tmp.png')
savefig('tmp.pdf')
raw_input()



