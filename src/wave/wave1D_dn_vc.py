#!/usr/bin/env python
"""
1D wave equation with Dirichlet or Neumann conditions
and variable wave velocity::

 u, x, t, cpu = solver(I, V, f, c, U_0, U_L, L, Nx, C, T,
                       user_action=None, version='scalar',
                       dt_safety_factor=1.0)

Solve the wave equation u_tt = (c**2*u_x)_x + f(x,t) on (0,L) with
u=U_0 or du/dn=0 on x=0, and u=u_L or du/dn=0
on x = L. If U_0 or U_L equals None, the du/dn=0 condition
is used, otherwise U_0(t) and/or U_L(t) are used for Dirichlet cond.
Initial conditions: u=I(x), u_t=V(x).

Nx is the total number of mesh intervals in x direction,
and mesh points are numbered from 0 to Nx.

C is the Courant number (=max(c)*dt/dx).
dt_safety_factor is used to computed the time step:
dt is the maximum value for stability times the dt_safety_factor.

T is the stop time for the simulation.

I, f, U_0, U_L, and c are functions: I(x), f(x,t), U_0(t),
U_L(t), c(x).
U_0 and U_L can also be 0, or None, where None implies
du/dn=0 boundary condition. f and V can also be 0 or None
(equivalent to 0). c can be a number or a function c(x).

user_action is a function of (u, x, t, n) where the calling code
can add visualization, error computations, data analysis,
store solutions, etc.
"""
import time, glob, shutil
from scitools.std import *

def solver(I, V, f, c, U_0, U_L, L, Nx, C, T,
           user_action=None, version='scalar',
           dt_safety_factor=1.0):
    """Solve u_tt=(c^2*u_x)_x + f on (0,L)x(0,T]."""
    import time;  t0 = time.clock()  # for measuring CPU time

    x = linspace(0, L, Nx+1)     # mesh points in space
    dx = x[1] - x[0]

    if isinstance(c, (float,int)):
        c = zeros(x.shape) + c
    else:
        c_ = zeros(x.shape)
        for i in range(Nx+1):
            c_[i] = c(x[i])
        c = c_

    dt = dt_safety_factor*C*dx/c.max()
    N = int(round(T/dt))
    t = linspace(0, N*dt, N+1)      # mesh points in time
    q = c**2
    C2 = (dt/dx)**2; dt2 = dt*dt    # help variables in the scheme

    # Wrap user-given f, V, U_0, U_L
    if f is None or f == 0:
        f = (lambda x, t: 0) if version == 'scalar' else \
            lambda x, t: zeros(x.shape)
    if V is None or V == 0:
        V = (lambda x: 0) if version == 'scalar' else \
            lambda x: zeros(x.shape)
    if U_0 is not None:
        if isinstance(U_0, (float,int)) and U_0 == 0:
            U_0 = lambda t: 0
    if U_L is not None:
        if isinstance(U_L, (float,int)) and U_L == 0:
            U_L = lambda t: 0

    u   = zeros(Nx+1)   # solution array at new time level
    u_1 = zeros(Nx+1)   # solution at 1 time level back
    u_2 = zeros(Nx+1)   # solution at 2 time levels back

    # Set initial condition
    for i in range(0,Nx+1):
        u_1[i] = I(x[i])

    if user_action is not None:
        user_action(u_1, x, t, 0)

    # Special formula for the first step
    for i in range(1, Nx):
        u[i] = u_1[i] + dt*V(x[i]) + \
        0.5*C2*(0.5*(q[i] + q[i+1])*(u_1[i+1] - u_1[i]) - \
                0.5*(q[i] + q[i-1])*(u_1[i] - u_1[i-1])) + \
        0.5*dt2*f(x[i], t[0])

    if U_0 is None:
        # Set boundary values (x=0: i-1 -> i+1 since u[i-1]=u[i+1]
        # when du/dn = 0, on x=L: i+1 -> i-1 since u[i+1]=u[i-1])
        i = 0
        ip1 = i+1
        im1 = ip1  # i-1 -> i+1
        u[i] = u_1[i] + dt*V(x[i]) + \
               0.5*C2*(0.5*(q[i] + q[ip1])*(u_1[ip1] - u_1[i])  - \
                       0.5*(q[i] + q[im1])*(u_1[i] - u_1[im1])) + \
        0.5*dt2*f(x[i], t[0])
    else:
        u[0] = U_0(dt)

    if U_L is None:
        i = Nx
        im1 = i-1
        ip1 = im1  # i+1 -> i-1
        u[i] = u_1[i] + dt*V(x[i]) + \
               0.5*C2*(0.5*(q[i] + q[ip1])*(u_1[ip1] - u_1[i])  - \
                       0.5*(q[i] + q[im1])*(u_1[i] - u_1[im1])) + \
        0.5*dt2*f(x[i], t[0])
    else:
        u[Nx] = U_L(dt)

    if user_action is not None:
        user_action(u, x, t, 1)

    u_2[:], u_1[:] = u_1, u

    for n in range(1, N):
        # Update all inner points
        if version == 'scalar':
            for i in range(1, Nx):
                u[i] = - u_2[i] + 2*u_1[i] + \
                    C2*(0.5*(q[i] + q[i+1])*(u_1[i+1] - u_1[i])  - \
                        0.5*(q[i] + q[i-1])*(u_1[i] - u_1[i-1])) + \
                dt2*f(x[i], t[n])

        elif version == 'vectorized':
            u[1:-1] = - u_2[1:-1] + 2*u_1[1:-1] + \
            C2*(0.5*(q[1:-1] + q[2:])*(u_1[2:] - u_1[1:-1]) -
                0.5*(q[1:-1] + q[:-2])*(u_1[1:-1] - u_1[:-2])) + \
            dt2*f(x[1:-1], t[n])
        else:
            raise ValueError('version=%s' % version)

        # Insert boundary conditions
        if U_0 is None:
            # Set boundary values
            # x=0: i-1 -> i+1 since u[i-1]=u[i+1] when du/dn=0
            # x=L: i+1 -> i-1 since u[i+1]=u[i-1] when du/dn=0
            i = 0
            ip1 = i+1
            im1 = ip1
            u[i] = - u_2[i] + 2*u_1[i] + \
                   C2*(0.5*(q[i] + q[ip1])*(u_1[ip1] - u_1[i])  - \
                       0.5*(q[i] + q[im1])*(u_1[i] - u_1[im1])) + \
            dt2*f(x[i], t[n])
        else:
            u[0] = U_0(t[n+1])

        if U_L is None:
            i = Nx
            im1 = i-1
            ip1 = im1
            u[i] = - u_2[i] + 2*u_1[i] + \
                   C2*(0.5*(q[i] + q[ip1])*(u_1[ip1] - u_1[i])  - \
                       0.5*(q[i] + q[im1])*(u_1[i] - u_1[im1])) + \
            dt2*f(x[i], t[n])
        else:
            u[Nx] = U_L(t[n+1])

        if user_action is not None:
            if user_action(u, x, t, n+1):
                break

        # Update data structures for next step
        u_2[:], u_1[:] = u_1, u

    cpu_time = t0 - time.clock()
    return u, x, t, cpu_time


import nose.tools as nt

def test_quadratic():
    """
    Check the scalar and vectorized versions work for
    a quadratic u(x,t)=x(L-x)(1+t) that is exactly reproduced,
    provided c(x) is constant.
    We simulate in [0, L/2] and apply a symmetry condition
    at the end x=L/2.
    """
    exact_solution = lambda x, t: x*(L-x)*(1+0.5*t)
    I = lambda x: exact_solution(x, 0)
    V = lambda x: 0.5*exact_solution(x, 0)
    f = lambda x, t: 2*(1+0.5*t)*c**2
    U_0 = lambda t: exact_solution(0, t)
    U_L = None
    L = 2.5
    c = 1.5
    Nx = 3  # very coarse mesh
    C = 1
    T = 18  # long time integration

    def assert_no_error(u, x, t, n):
        u_e = exact_solution(x, t[n])
        diff = abs(u - u_e).max()
        #print n, diff
        #print u
        #print u_e
        nt.assert_almost_equal(diff, 0, places=13)

    solver(I, V, f, c, U_0, U_L, L/2, Nx, C, T,
           user_action=assert_no_error, version='scalar',
           dt_safety_factor=1)
    solver(I, V, f, c, U_0, U_L, L/2, Nx, C, T,
           user_action=assert_no_error, version='vectorized',
           dt_safety_factor=1)

def test_plug():
    """Check that an initial plug is correct back after one period."""
    L = 1
    I = lambda x: 0 if abs(x-L/2.0) > 0.1 else 1

    u_s, x, t, cpu = solver(
        I=I,
        V=None, f=None, c=0.5, U_0=None, U_L=None, L=L,
        Nx=50, C=1, T=4, user_action=None, version='scalar')
    u_v, x, t, cpu = solver(
        I=I,
        V=None, f=None, c=0.5, U_0=None, U_L=None, L=L,
        Nx=50, C=1, T=4, user_action=None, version='vectorized')
    diff = abs(u_s - u_v).max()
    nt.assert_almost_equal(diff, 0, places=13)
    u_0 = array([I(x_) for x_ in x])
    diff = abs(u_s - u_0).max()
    nt.assert_almost_equal(diff, 0, places=13)


class PlotSolution:
    """
    Class for the user_action function in solver.
    Visualizes the solution only.
    """
    def __init__(self,
                 casename='tmp',    # prefix in filenames
                 umin=-1, umax=1,   # fixed range of y axis
                 pause_between_frames=None,  # movie speed
                 backend='matplotlib',       # or 'gnuplot'
                 screen_movie=True, # show movie on screen?
                 every_frame=1):    # show every_frame frame
        self.casename = casename
        self.yaxis = [umin, umax]
        self.pause = pause_between_frames
        module = 'scitools.easyviz.' + backend + '_'
        exec('import %s as st' % module)
        self.st = st
        self.screen_movie = screen_movie
        self.every_frame = every_frame

        # Clean up old movie frames
        for filename in glob('frame_*.png'):
            os.remove(filename)

    def __call__(self, u, x, t, n):
        if n % self.every_frame != 0:
            return
        self.st.plot(x, u, 'r-',
                     xlabel='x', ylabel='u',
                     axis=[x[0], x[-1],
                           self.yaxis[0], self.yaxis[1]],
                     title='t=%f' % t[n],
                     show=self.screen_movie)
        # pause
        if t[n] == 0:
            time.sleep(2)  # let initial condition stay 2 s
        else:
            if self.pause is None:
                pause = 0.2 if u.size < 100 else 0
            time.sleep(pause)

        self.st.savefig('%s_frame_%04d.png' % (self.casename, n))

    def make_movie_file(self):
        """
        Create subdirectory based on casename, move all plot
        frame files to this directory, and generate
        an index.html for viewing the movie in a browser
        (as a sequence of PNG files).
        """
        directory = self.casename
        shutil.rmtree(directory)   # rm -rf directory
        os.mkdir(directory)        # mkdir directory
        # mv frame_*.png directory
        for filename in glob('frame_*.png'):
            os.rename(filename, os.path.join(directory, filename))
        os.chdir(directory)        # cd directory
        self.st.movie('frame_*.png', encoder='html',
                      output_file='index.html', fps=4)
        # Can make other movie files...

def moving_end(C=1, Nx=50, reflecting_right_boundary=True,
               version='vectorized'):
    L = 1.
    c = 1
    T = 2
    I = lambda x: 0

    def U_0(t):
        return 0.25*sin(6*pi*t) if t < 0.75 else 0

    if reflecting_right_boundary:
        U_L = None
    else:
        U_L = 0

    action = PlotSolution('moving_end', -1, 1)
    solver(I, V, f, c, U_0, U_L, L, Nx, C, T,
           user_action=action, version=version,
           dt_safety_factor=1)


class PlotMediumAndSolution(PlotSolution):
    def __init__(self, medium, **kwargs):
        """Mark medium in plot: medium=[x_L, x_R]."""
        self.medium = medium
        PlotSolution.__init__(self, **kwargs)

    def __call__(self, u, x, t, n):
        # Plot u and mark medium x=x_L and x=x_R
        x_L, x_R = self.medium
        umin, umax = self.yaxis
        self.st.plot(x, u, 'r-',
                     [x_L, x_L], [umin, umax], 'k--',
                     [x_R, x_R], [umin, umax], 'k--',
                     xlabel='x', ylabel='u',
                     axis=[x[0], x[-1], umin, umax],
                     title='Nx=%d, t=%f' % (x.size-1, t[n]))
        if t[n] == 0:
            time.sleep(2)  # let initial condition stay 2 s
        # No sleep - this is used for large meshes
        self.st.savefig('frame_%04d.png' % n)  # for movie making


def pulse(C=1, Nx=200, animate=True, version='vectorized', T=2,
          loc='center', pulse_tp='gaussian', slowness_factor=2,
          medium=[0.7, 0.9], every_frame=1):
    """
    Various peaked-shaped initial conditions on [0,1].
    Wave velocity is decreased by the slowness_factor inside
    the medium spcification. The loc parameter can be 'center'
    or 'left', depending on where the pulse is to be located.
    """
    L = 1.
    if loc == 'center':
        xc = L/2
    elif loc == 'left':
        xc = 0
    sigma = L/20.  # width measure of I(x)
    c_0 = 1.0      # wave velocity outside medium

    if pulse_tp in ('gaussian','Gaussian'):
        def I(x):
            return exp(-0.5*((x-xc)/sigma)**2)
    elif pulse_tp == 'plug':
        def I(x):
            return 0 if abs(x-xc) > sigma else 1
    elif pulse_tp == 'cosinehat':
        def I(x):
            # One period of a cosine
            w = 2
            a = w*sigma
            return 0.5*(1 + cos(pi*(x-xc)/a)) \
                   if xc - a <= x <= xc + a else 0

    elif pulse_tp == 'half-cosinehat':
        def I(x):
            # Half a period of a cosine
            w = 4
            a = w*sigma
            return cos(pi*(x-xc)/a) \
                   if xc - 0.5*a <= x <= xc + 0.5*a else 0
    else:
        raise ValueError('Wrong pulse_tp="%s"' % pulse_tp)

    def c(x):
        return c_0/slowness_factor \
               if medium[0] <= x <= medium[1] else c_0

    umin=-0.5; umax=1.5*I(xc)
    casename = '%s_Nx%s_sf%s' % \
               (pulse_tp, Nx, slowness_factor)
    action = PlotMediumAndSolution(
        medium, casename=casename, umin=umin, umax=umax,
        every_frame=every_frame, screen_movie=animate)

    solver(I=I, V=None, f=None, c=c, U_0=None, U_L=None,
           L=L, Nx=Nx, C=C, T=T,
           user_action=action, version=version,
           dt_safety_factor=1)



if __name__ == '__main__':
    import sys
    from scitools.misc import function_UI
    cmd = function_UI([test_quadratic, test_plug, pulse,
                       moving_end,], sys.argv)
    eval(cmd)
    raw_input()
