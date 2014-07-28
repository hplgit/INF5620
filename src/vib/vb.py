from numpy import *
from matplotlib.pyplot import *

def solver(I, V, m, c, s, F, dt, T):
    """
    Solve m*u'' + cu' + s(u) = F(t) for t in (0,T],
    u(0)=I and u'(0)=V,
    by a central finite difference method with time step dt.
    """
    dt = float(dt)
    N = int(round(T/dt))
    u = zeros(N+1)
    t = linspace(0, N*dt, N+1)

    u[0] = I
    u[1] = u[0] + dt*(c/m*dt - 1)*V + \
           1./(2*m)*dt**2*(F(t[0]) - s(u[0]))
    for n in range(1,N):
        u[n+1] = 1./(m + c*dt)*(2*m*u[n] + \
                 (c*dt - m)*u[n-1] + dt**2/m*(F(t[n]) - s(u[n])))
    return u, t

def visualize(u, t):
    plot(t, u, 'b-')
    xlabel('t')
    ylabel('u')
    dt = t[1] - t[0]
    title('dt=%g' % dt)
    umin = 1.2*u.min(); umax = -umin
    axis([t[0], t[-1], umin, umax])
    savefig('vib2.png')
    savefig('vib2.pdf')
    savefig('vib2.eps')

import nose.tools as nt

def test_constant_constant():
    """Verify a constant solution."""

    def exact_solution(t):
        return I

    I = 1.2; V = 0; m = 2; c = 0.9
    w = 1.5
    s = lambda u: w**2*u
    F = lambda t: w**2*exact_solution(t)
    dt = 0.2
    T = 2
    u, t = solver(I, V, m, c, s, F, dt, T)
    difference = abs(exact_solution(t) - u).max()
    nt.assert_almost_equal(difference, 0, places=14)

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--I', type=float, default=1.0)
    parser.add_argument('--V', type=float, default=0.0)
    parser.add_argument('--m', type=float, default=1.0)
    parser.add_argument('--c', type=float, default=0.0)
    parser.add_argument('--s', type=str, default='u')
    parser.add_argument('--F', type=str, default='0')
    parser.add_argument('--dt', type=float, default=0.05)
    parser.add_argument('--T', type=float, default=140)
    parser.add_argument('--window_width', type=float, default=30)
    parser.add_argument('--savefig', action='store_true')
    # Hack to allow --SCITOOLS options (read when importing scitools.std)
    parser.add_argument('--SCITOOLS_easyviz_backend', default='matplotlib')
    a = parser.parse_args()
    from scitools.std import StringFunction
    s = StringFunction(a.s, independent_variable='u')
    F = StringFunction(a.F, independent_variable='t')
    I, V, m, c, dt, T, window_width, savefig = \
       a.I, a.V, a.m, a.c, a.dt, a.T, a.window_width, a.savefig

    u, t = solver(I, V, m, c, s, F, dt, T)
    num_periods = empirical_freq_and_amplitude(u, t)
    if num_periods <= 40:
        figure()  # new matplotlib figure
        visualize(u, t)
    else:
        visualize_front(u, t, window_width, savefig)  # scitools
        visualize_front_ascii(u, t)
    show()

def empirical_freq_and_amplitude(u, t):
    minima, maxima = minmax(t, u)
    p = periods(maxima)
    a = amplitudes(minima, maxima)
    figure()
    w = 2*pi/p
    plot(range(len(p)), w, 'r-')
    hold('on')
    plot(range(len(a)), a, 'b-')
    ymax = 1.1*max(w.max(), a.max())
    ymin = 0.9*min(w.min(), a.min())
    axis([0, max(len(p), len(a)), ymin, ymax])
    legend(['estimated frequency', 'estimated amplitude'],
           loc='upper right')
    return len(maxima)

def visualize_front(u, t, window_width, savefig=False):
    """
    Visualize u and the exact solution vs t, using a
    moving plot window and continuous drawing of the
    curves as they evolve in time.
    Makes it easy to plot very long time series.
    P is the approximate duration of one period.
    """
    import scitools.std as st
    from scitools.MovingPlotWindow import MovingPlotWindow

    umin = 1.2*u.min();  umax = -umin
    plot_manager = MovingPlotWindow(
        window_width=window_width,
        dt=t[1]-t[0],
        yaxis=[umin, umax],
        mode='continuous drawing')
    for n in range(1,len(u)):
        if plot_manager.plot(n):
            s = plot_manager.first_index_in_plot
            st.plot(t[s:n+1], u[s:n+1], 'r-1',
                    title='t=%6.3f' % t[n],
                    axis=plot_manager.axis(),
                    show=not savefig) # drop window if savefig
            if savefig:
                print 't=%g' % t[n]
                st.savefig('tmp_vib%04d.png' % n)
        plot_manager.update(n)

def visualize_front_ascii(u, t, fps=10):
    """
    Plot u and the exact solution vs t line by line in a
    terminal window (only using ascii characters).
    Makes it easy to plot very long time series.
    """
    from scitools.avplotter import Plotter
    import time
    umin = 1.2*u.min();  umax = -umin

    p = Plotter(ymin=umin, ymax=umax, width=60, symbols='+o')
    for n in range(len(u)):
        print p.plot(t[n], u[n]), '%.2f' % (t[n])
        time.sleep(1/float(fps))

def minmax(t, u):
    """
    Compute all local minima and maxima of the function u(t),
    represented by discrete points in the arrays u and t.
    Return lists minima and maxima of (t[i],u[i]) extreme points.
    """
    minima = []; maxima = []
    for n in range(1, len(u)-1, 1):
        if u[n-1] > u[n] < u[n+1]:
            minima.append((t[n], u[n]))
        if u[n-1] < u[n] > u[n+1]:
            maxima.append((t[n], u[n]))
    return minima, maxima

def periods(extrema):
    """
    Given a list of (t,u) points of the maxima or minima,
    return an array of the corresponding local periods.
    """
    p = [extrema[n][0] - extrema[n-1][0]
         for n in range(1, len(extrema))]
    return array(p)

def amplitudes(minima, maxima):
    """
    Given a list of (t,u) points of the minima and maxima of
    u, return an array of the corresponding local amplitudes.
    """
    # Compare first maxima with first minima and so on
    a = [(abs(maxima[n][1] - minima[n][1]))/2.0
         for n in range(min(len(minima),len(maxima)))]
    return array(a)

if __name__ == '__main__':
    main()
    #test_constant()

