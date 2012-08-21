from numpy import *

def theta_rule(I, a, T, dt, theta):
    """Solve u'=-a*u, u(0)=I, for t in (0,T] with steps of dt."""
    N = int(round(T/float(dt)))  # no of intervals
    u = zeros(N+1)
    t = linspace(0, T, N+1)

    u[0] = I
    for n in range(0, N):
        u[n+1] = (1 - (1-theta)*a*dt)/(1 + theta*dt*a)*u[n]
    return u, t

def verify_three_steps():
    # Three manual steps
    theta = 0.8; a = 2; I = 0.1; dt = 0.8
    factor = (1 - (1-theta)*a*dt)/(1 + theta*dt*a)
    u1 = factor*I
    u2 = factor*u1
    u3 = factor*u2

    N = 3  # number of time steps
    u, t = theta_rule(I=I, a=a, T=N*dt, dt=dt, theta=theta)

    print u[1:]  # u[1], u[2], ...
    print u1, u2, u3

# Better version:

def verify_three_steps():
    # Three manual steps
    theta = 0.8; a = 2; I = 0.1; dt = 0.8
    factor = (1 - (1-theta)*a*dt)/(1 + theta*dt*a)
    u1 = factor*I
    u2 = factor*u1
    u3 = factor*u2

    N = 3  # number of time steps
    u, t = theta_rule(I=I, a=a, T=N*dt, dt=dt, theta=theta)

    tol = 1E-15  # tolerance for comparing floats
    difference = abs(u1-u[1]) + abs(u2-u[2]) + abs(u3-u[3])
    success = difference <= tol
    return success

def main():
    u, t = theta_rule(I=1, a=2, T=8, dt=0.8, theta=1)
    # Write out a table of t and u values:
    for i in range(len(t)):
        print 't=%6.3f u=%g' % (t[i], u[i])
        # or print 't={t:6.3f} u={u:g}'.format(t=t[i], u=u[i])

if verify_three_steps():
    main()
else:
    print 'Bug in the implementation!'




