from scitools.std import *

def theta_rule(I, a, T, dt, theta):
    """Solve u'=-a*u, u(0)=I, for t in (0,T] with steps of dt."""
    N = int(round(T/float(dt)))  # no of intervals
    u = zeros(N+1)
    t = linspace(0, T, N+1)

    u[0] = I
    for n in range(0, N):
        u[n+1] = (1 - (1-theta)*a*dt)/(1 + theta*dt*a)*u[n]
    return u, t


def exact_solution(t, I, a):
    return I*exp(-a*t)

def explore(I, a, T, dt, theta=0.5, makeplot=True):
    """
    Run a case with the theta_rule, compute error measure,
    and plot the numerical and exact solutions (if makeplot=True).
    """
    u, t = theta_rule(I, a, T, dt, theta)  # Numerical solution
    u_e = exact_solution(t, I, a)
    e = u_e - u
    E = sqrt(dt*sum(e**2))
    if makeplot:
        figure()                         # create new plot
        t_e = linspace(0, T, 1001)       # very fine mesh for u_e
        u_e = exact_solution(t_e, I, a)
        theta2name = {0: 'FE', 1: 'BE', 0.5: 'CN'}
        plot(t,   u,   'ro',             # red circles for u
             t_e, u_e, 'b-',             # blue line for u_e
             legend=['numerical', 'exact'],
             xlabel='t',
             ylabel='u',
             title='Method: theta-rule, theta=%g, dt=%g' %
             (theta, dt),
             savefig='%s_%g.png' % (theta2name[theta], dt),
             show=True)
    return E

I = 1
a = 2
T = 5
for theta in 0, 0.5, 1:
    for dt in 0.4, 0.04:
        E = explore(I, a, T, dt, theta, makeplot=True)
        print '%3.1f %6.2f: %12.3E' % (theta, dt, E)
