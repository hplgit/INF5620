
from scitools.MovingPlotWindow import MovingPlotWindow
import time


def central_nondamped(I, w, dt, T):
    """
    Solve vibrating system equation u'' + w**2 = 0 for t
    in (0,T], u(0)=I and u'(0)=0, by a central finite
    difference method with time step dt. Plot the
    solution and return u, t.
    """

    if dt > 2./w:
        print 'Unstable scheme'
    N = int(round(T/float(dt)))
    u = zeros(N+1)
    t = linspace(0, T, N+1)

    umin = -1.2*I
    umax = -umin
    period = 2*pi/w
    plot_manager = MovingPlotWindow(
        window_width=8*period,
        dt=dt,
        yaxis=[umin, umax],
        mode='continuous drawing')
    u[0] = I
    u[1] = u[0] - 0.5*dt**2*w**2*u[0]
    for n in range(1,N):
        u[n+1] = 2*u[n] - u[n-1] - dt**2*w**2*u[n]

        if plot_manager.plot(n):
            s = plot_manager.first_index_in_plot
            plot(t[s:n+2], u[s:n+2], 'r-',
                 t[s:n+2], I*cos(w*t)[s:n+2], 'b-',
                 title='t=%6.3f' % t[n+1],
                 axis=plot_manager.axis())
        plot_manager.update(n)
    return u, t

def central_damped(I, w, b, F, dt, T):
    """
    Solve a vibrating system problem u'' + b*u' + w**2*u = F(t)
    with u(0)=I and u'(0)=0, for t in (0,T], by a central
    finite difference method with time step dt.
    Plot the solution with a moving window, and return u, t.
    """
    if dt > 2./w:
        print 'Unstable scheme'
    N = int(round(T/float(dt)))
    u = zeros(N+1)
    t = linspace(0, T, N+1)

    umin = -1.2*I
    umax = -umin
    period = 2*pi/w
    plot_manager = MovingPlotWindow(8*period, dt, yaxis=[umin, umax],
                                    mode='continuous drawing')

    Cp = b/2*dt + 1
    Cm = b/2*dt - 1
    D = 2 - dt**2*w**2

    u[0] = I
    u[1] = 0.5*D*u[0] + 0.5*dt**2*F(0)

    for n in range(1,N):
        u[n+1] = (Cm*u[n-1] + D/Cp*u[n] + dt**2*F(n*dt))/Cp

        if plot_manager.plot(n):
            s = plot_manager.first_index_in_plot
            plot(t[s:n+2], u[s:n+2], 'r-',
                 title='t=%6.3f' % t[n+1],
                 axis=plot_manager.axis())
        plot_manager.update(n)
    return u, t

#central_nondamped(1, pi, 0.25, 400)

#plot_c()

central_damped(I=1, w=1, b=0.05, F=lambda t: 2*sin(4*t), dt=0.05, T=20)















# calc damped
# write analysis undamped
# write analysis damped
# write derivation
# program full dumpy road
# figure for bumpy road
# 1D wave eq u=UL(t),UR(t) BC
# same du/dn=0 BC
# 2D wave eq w/du/dn=0 and plotting
# 1D wave eq in spherical symmetry
# solve ODE in SymPy, sympy.dsolve
# ODE.py, liwei inspired ODESolver with all exerc from TCSE6 + new schemes


