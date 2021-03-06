#!/usr/bin/env python
"""
1D wave equation with u=0 at the boundary.
The solver function here offers scalar and vectorized versions.
See wave1D_u0_s.py for documentation. The only difference
is that function solver takes an additional argument "version":
version='scalar' implies explicit loops over mesh point,
while version='vectorized' provides a vectorized version.
"""
from numpy import *

def solver(I, V, f, c, L, Nx, C, T, user_action=None,
           version='vectorized'):
    """Solve u_tt=c^2*u_xx + f on (0,L)x(0,T]."""
    x = linspace(0, L, Nx+1)   # mesh points in space
    dx = x[1] - x[0]
    dt = C*dx/c
    N = int(round(T/dt))
    t = linspace(0, N*dt, N+1) # mesh points in time
    C2 = C**2                  # help variable in the scheme
    if f is None or f == 0:
        f = lambda x, t: 0 if version == 'scalar' else \
            lambda x, t: zeros(x.shape)
    if V is None or V == 0:
        V = lambda x: 0 if version == 'scalar' else \
            lambda x: zeros(x.shape)

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

        if version == 'scalar':
            for i in range(1, Nx):
                u[i] = - u_2[i] + 2*u_1[i] + \
                       C2*(u_1[i-1] - 2*u_1[i] + u_1[i+1]) + \
                       dt**2*f(x[i], t[n])
        elif version == 'vectorized':
            f_a = f(x, t[n])  # precompute in array
            u[1:-1] = - u_2[1:-1] + 2*u_1[1:-1] + \
                C2*(u_1[0:-2] - 2*u_1[1:-1] + u_1[2:]) + \
                dt**2*f_a[1:-1]
        elif version == 'vectorized2':
            f_a = f(x, t[n])  # precompute in array
            u[1:Nx] =  - u_2[1:Nx] + 2*u_1[1:Nx] + \
                C2*(u_1[0:Nx-1] - 2*u_1[1:Nx] + u_1[2:Nx+1]) + \
                dt**2*f_a[1:Nx]

        # Insert boundary conditions
        u[0] = 0;  u[Nx] = 0
        if user_action is not None:
            if user_action(u, x, t, n+1):
                break

        # Switch variables before next step
        u_2[:], u_1[:] = u_1, u

    cpu_time = t0 - time.clock()
    return u, x, t, cpu_time

def viz(I, V, f, c, L, Nx, C, T, umin, umax, animate=True,
        version='vectorized'):
    """Run solver and visualize u at each time level."""
    import scitools.std as st, time, glob, os
    #num_frames = 100 # max no of frames in movie

    def plot_u(u, x, t, n):
        """user_action function for solver."""
        try:
            every = t.size/num_frames
        except NameError:
            every = 1  # plot every frame
        if n % every == 0:
            st.plot(x, u, 'r-',
                    xlabel='x', ylabel='u',
                    axis=[0, L, umin, umax],
                    title='t=%f' % t[n])
            # Let the initial condition stay on the screen for 2
            # seconds, else insert a pause of 0.2 s between each plot
            time.sleep(2) if t[n] == 0 else time.sleep(0.2)
            st.savefig('frame_%04d.png' % n)  # for movie making

    # Clean up old movie frames
    for filename in glob.glob('frame_*.png'):
        os.remove(filename)

    user_action = plot_u if animate else None
    u, x, t, cpu = solver(I, V, f, c, L, Nx, C, T,
                          user_action, version)
    # Make movie files
    st.movie('frame_*.png', encoder='mencoder', fps=4,
             output_file='movie.avi')
    st.movie('frame_*.png', encoder='html', fps=4,
             output_file='movie.html')
    return cpu

import nose.tools as nt

def test_quadratic():
    """
    Check the scalar and vectorized versions work for
    a quadratic u(x,t)=x(L-x)(1+t) that is exactly reproduced.
    """
    exact_solution = lambda x, t: x*(L - x)*(1 + 0.5*t)
    I = lambda x: exact_solution(x, 0)
    V = lambda x: 0.5*exact_solution(x, 0)
    f = lambda x, t: 2*c**2*(1 + 0.5*t)
    L = 2.5
    c = 1.5
    Nx = 3  # very coarse mesh
    C = 1
    T = 18  # long time integration

    def assert_no_error(u, x, t, n):
        u_e = exact_solution(x, t[n])
        diff = abs(u - u_e).max()
        print diff
        nt.assert_almost_equal(diff, 0, places=13)

    solver(I, V, f, c, L, Nx, C, T,
           user_action=assert_no_error, version='scalar')
    solver(I, V, f, c, L, Nx, C, T,
           user_action=assert_no_error, version='vectorized')

def guitar():
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
    C = 0.85

    def I(x):
        return a*x/x0 if x < x0 else a/(L-x0)*(L-x)

    umin = -1.2*a;  umax = -umin
    cpu = viz(I, 0, 0, c, L, Nx, C, T, umin, umax, animate=True)

def run_efficiency_experiments():
    L = 1
    x0 = 0.8*L
    a = 1
    c = 2
    T = 8
    C = 0.9
    umin = -1.2*a;  umax = -umin

    def I(x):
        return a*x/x0 if x < x0 else a/(L-x0)*(L-x)

    for Nx in [50, 100, 200, 400, 800]:
        print 'solving scalar Nx=%d' % Nx,
        cpu_s = viz(I, 0, 0, c, L, Nx, C, T, umin, umax,
                    animate=False, version='scalar')
        print cpu_s
        print 'solving vectorized Nx=%d' % Nx,
        cpu_v = viz(I, 0, 0, c, L, Nx, C, T, umin, umax,
                    animate=False, version='vectorized')
        print cpu_v
        print 'Nx=%3d: cpu_v/cpu_s: %.3f' % (Nx, cpu_v/float(cpu_s))

if __name__ == '__main__':
    test_quadratic()  # verify
    guitar()
    #run_efficiency_experiments()

