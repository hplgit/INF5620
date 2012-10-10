#!/usr/bin/env python
"""
1D wave equation with u=0 at the boundary.
Simplest possible implementation.

The key function is::

  u, x, t, cpu = (I, V, f, c, L, Nx, C, T, user_action)

which solves the wave equation u_tt = c**2*u_xx on (0,L) with u=0
on x=0,L, for t in (0,T].  Initial conditions: u=I(x), u_t=V(x).

Nx is the total number of spatial mesh cells; mesh points
are numbered from 0 to Nx.
C is the Courant number (=c*dt/dx), which specifies dt.
T is the stop time for the simulation.
f(x,t) is a function for the source term (can be 0 or None).
I and V are functions of x.

user_action is a function of (u, x, t, n) where the calling
code can add visualization, error computations, etc.
"""

from numpy import *

def solver(I, V, f, c, L, Nx, C, T, user_action=None):
    """Solve u_tt=c^2*u_xx + f on (0,L)x(0,T]."""
    x = linspace(0, L, Nx+1)   # mesh points in space
    dx = x[1] - x[0]
    dt = C*dx/c
    N = int(round(T/dt))
    t = linspace(0, N*dt, N+1) # mesh points in time
    C2 = C**2                  # help variable in the scheme
    if f is None or f == 0 :
        f = lambda x, t: 0
    if V is None or V == 0:
        V = lambda x: 0

    u   = zeros(Nx+1)   # solution array at new time level
    u_1 = zeros(Nx+1)   # solution at 1 time level back
    u_2 = zeros(Nx+1)   # solution at 2 time levels back

    import time;  t0 = time.clock()  # for measuring CPU time

    # Load initial condition into u_1
    for i in range(0,Nx+1):
        u_1[i] = I(x[i])

    if user_action is not None:
        user_action(u_1, x, t, 0)

    # Special formula for first time step
    n = 0
    for i in range(1, Nx):
        u[i] = u_1[i] + dt*V(x[i]) + \
               0.5*C2*(u_1[i-1] - 2*u_1[i] + u_1[i+1]) + \
               0.5*dt**2*f(x[i], t[n])
    u[0] = 0;  u[Nx] = 0

    if user_action is not None:
        user_action(u, x, t, 1)

    u_2[:], u_1[:] = u_1, u

    for n in range(1, N):
        # Update all inner points at time t[n+1]
        for i in range(1, Nx):
            u[i] = - u_2[i] + 2*u_1[i] + \
                     C2*(u_1[i-1] - 2*u_1[i] + u_1[i+1]) + \
                     dt**2*f(x[i], t[n])

        # Insert boundary conditions
        u[0] = 0;  u[Nx] = 0
        if user_action is not None:
            if user_action(u, x, t, n+1):
                break

        # Switch variables before next step
        u_2[:], u_1[:] = u_1, u

    cpu_time = t0 - time.clock()
    return u, x, t, cpu_time

import nose.tools as nt

def test_quadratic():
    """Check that u(x,t)=x(L-x)(1+t) is exactly reproduced."""
    def exact_solution(x, t):
        return x*(L-x)*(1 + 0.5*t)

    def I(x):
        return exact_solution(x, 0)

    def V(x):
        return 0.5*exact_solution(x, 0)

    def f(x, t):
        return 2*(1 + 0.5*t)*c**2

    L = 2.5
    c = 1.5
    Nx = 3  # very coarse mesh
    C = 0.75
    T = 18

    u, x, t, cpu = solver(I, V, f, c, L, Nx, C, T)
    u_e = exact_solution(x, t[-1])
    diff = abs(u - u_e).max()
    nt.assert_almost_equal(diff, 0, places=14)

def viz(I, V, f, c, L, Nx, C, T, umin, umax, animate=True):
    """Run solver and visualize u at each time level."""
    import scitools.std as st, time, glob, os

    def plot_u(u, x, t, n):
        """user_action function for solver."""
        st.plot(x, u, 'r-',
                xlabel='x', ylabel='u',
                axis=[0, L, umin, umax],
                title='t=%f' % t[n], show=True)
        # Let the initial condition stay on the screen for 2
        # seconds, else insert a pause of 0.2 s between each plot
        time.sleep(2) if t[n] == 0 else time.sleep(0.2)
        st.savefig('frame_%04d.png' % n)  # for movie making

    # Clean up old movie frames
    for filename in glob.glob('frame_*.png'):
        os.remove(filename)

    user_action = plot_u if animate else None
    u, x, t, cpu = solver(I, V, f, c, L, Nx, C, T, user_action)

    # Make movie files
    st.movie('frame_*.png', encoder='mencoder', fps=4,
             output_file='movie.avi')
    st.movie('frame_*.png', encoder='html', fps=4,
             output_file='movie.html')

def guitar(C):
    """Triangular wave (pulled guitar string)."""
    L = 0.4
    x0 = 0.8*L
    a = 0.005
    freq = 440
    wavelength = 2*L
    c = freq*wavelength
    omega = 2*pi*freq
    num_periods = 1
    T = 2*pi/omega*num_periods
    Nx = 50

    def I(x):
        return a*x/x0 if x < x0 else a/(L-x0)*(L-x)

    umin = -1.2*a;  umax = -umin
    cpu = viz(I, 0, 0, c, L, Nx, C, T, umin, umax, animate=True)


if __name__ == '__main__':
    test_quadratic()
    import sys
    try:
        C = float(sys.argv[1])
        print 'C=%g' % C
    except IndexError:
        C = 0.85
    guitar(C)
