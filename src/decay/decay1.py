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

u, t = theta_rule(I=1, a=2, T=8, dt=0.8, theta=1)
# Write out a table of t and u values:
for i in range(len(t)):
    print 't=%6.3f u=%g' % (t[i], u[i])
    # or print 't={t:6.3f} u={u:g}'.format(t=t[i], u=u[i])

