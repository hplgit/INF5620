import numpy as np
import matplotlib.pyplot as plt
import sys, os
import nose.tools as nt

def leapfrog(I, a, b, T, dt, first_step='ForwardEuler'):
    dt = float(dt)
    Nt = int(round(T/dt))  # no of intervals
    u = np.zeros(Nt+1)
    t = np.linspace(0, Nt*dt, Nt+1)

    n = 0
    u[n] = I
    if first_step == 'ForwardEuler':
        u[n+1] = (1 - dt*a(t[n]))*u[n] + dt*b(t[n])
    elif first_step == 'CrankNicolson':
        u[n+1] = ((1 - 0.5*a(t[n])*dt)*u[n] + \
                  0.5*dt*(b(t[n]) + b(t[n+1])))/\
                  (1 + 0.5*dt*a(t[n+1]))
    else:
        raise ValueError('Wrong first_step=%s' % first_step)

    for n in range(1, Nt):
        u[n+1] = u[n-1] - 2*dt*a(t[n])*u[n] + 2*dt*b(t[n])
    return u, t

def visualize(u, t, t_e=None, u_e=None, show=True, title=None):
    """Plot u and (if available) the exact solution u_e, t_e."""
    plt.figure()                                # create new plot
    plt.plot(t, u, 'r-', label='numerical')     # red line
    if u_e is not None and t_e is not None:
        plt.plot(t_e, u_e, 'b-', label='exact') # blue line for u_e
    plt.xlabel('t')
    plt.ylabel('u')
    if title is None:
        plt.title('Leapfrog scheme')
    else:
        plt.title(title)
    plt.savefig('tmp_leapfrog.png')
    if show:
        plt.show()

def first_test_plot():
    u, t = leapfrog(I=1, a=lambda t: 1, b=lambda t: 0, T=4,
                    dt=0.01, first_step='ForwardEuler')
    visualize(u, t, u_e=np.exp(-t), t_e=t)

import sympy as sm

def analyze(u):
    """For given u(t) (Python function of t), perform MMS and
    show residual in the discrete equations."""

    t, dt, a = sm.symbols('t dt a')

    print 'Analyzing u_e(t)=%s' % u(t)
    print 'u(0)=%s' % u(t).subs(t, 0)

    # Fit source term to the given u(t)
    b = sm.diff(u(t), t) + a*u(t)
    b = sm.simplify(b)
    print 'Source term b:', b

    # Residual in discrete equations; Forward Euler step
    R_step1 = (u(t+dt) - u(t))/dt + a*u(t) - b
    R_step1 = sm.simplify(R_step1)
    print 'Residual Forward Euler step:', R_step1

    # Residual in discrete equations; Crank-Nicolson step
    R_step1 = (u(t+dt) - u(t))/dt + sm.Rational(1,2)*(
        a*u(t) + a*u(t+dt) - (b + b.subs(t, t+dt)))
    R_step1 = sm.simplify(R_step1)
    print 'Residual Crank-Nicolson step:', R_step1

    # Residual in discrete equations; Leapfrog steps
    R = (u(t+dt) - u(t-dt))/(2*dt) + a*u(t) - b
    R = sm.simplify(R)
    print 'Residual Leapfrog steps:', R


def test_linear_solution():
    """Nose test of exact reproduction of a linear solution."""
    c = 1.2
    I = 4.2

    def u_e(t):
        return c*t + I

    a = lambda t: 0.2  # just a constant
    b = lambda t: a(t)*u_e(t) + c

    T = 20
    #T = 200  # instability
    dt = 0.5
    u, t = leapfrog(I, a, b, T, dt, first_step='ForwardEuler')
    #visualize(u, t)
    error = u_e(t) - u
    max_error = np.abs(error).max()
    nt.assert_almost_equal(max_error, 0, delta=1E-13)

def convergence_rates():
    """Compute convergence rates for a given manufactured solution."""
    u_e = lambda t: np.sin(t)
    a = lambda t: 0.2
    b = lambda t: a(t)*np.sin(t) + np.cos(t)

    I = u_e(0)
    T = 4*np.pi
    dt_value = 0.5
    dt = []
    E = []
    for num_meshes in range(6):
        u, t = leapfrog(I, a, b, T, dt_value, first_step='ForwardEuler')
        #u, t = leapfrog(I, a, b, T, dt_value, first_step='CrankNicolson')
        #visualize(u, t)
        error_l2 = np.sqrt(dt_value*np.sum((u_e(t) - u)**2))
        #error_l2 = np.abs(u_e(t) - u).max()
        E.append(error_l2)
        dt.append(dt_value)
        dt_value /= 2.      # half dt for each mesh

    print 'Convergence rates:'
    print 'E:',  E
    print 'dt:', dt
    r = [np.log(E[i]/E[i-1])/np.log(dt[i]/dt[i-1])
         for i in range(1, len(E))]
    print 'r:', r
    return r

def test_convergence_rate():
    r = convergence_rates()
    nt.assert_almost_equal(r[-1], 2.0, places=1)

def demonstrate_instability():
    """Run larger and larger dt and demonstrate small, growing instability."""
    a = lambda t: 1.0
    b = lambda t: 0
    I = 1.0
    # Solution: exp(-t)
    T = 8
    dt = 0.001
    for k in range(8):
        u, t = leapfrog(I, a, b, T, dt, first_step='ForwardEuler')
        visualize(u, t, title='timestep: %g' % dt)
        dt *= 2

    # Even with very small steps we get instability after sufficiently
    # long time
    dt = 1E-3
    T = 15
    u, t = leapfrog(I, a, b, T, dt, first_step='ForwardEuler')
    visualize(u, t, title='timestep: %g' % dt)


def exact_discrete_solution():
    """Find exact discrete (analytical) solution. Use this
    to explain the observed instability."""
    A, n, p = sm.symbols('A n p')

    # u[n+1] = u[n-1] - 2*dt*a*u[n]
    # p = a*dt
    # u[n+1] - u[n-1] + 2*p*u[n] = 0
    # u[n] = A**n
    # Insert this u:
    R = A**(n+1) - A**(n-1) + 2*p*A**n
    R = R/A**(n-1)
    R = sm.simplify(R)
    print 'polynomial equation in A:', R, '= 0'

    print
    # Solve equation for A
    A1, A2 = sm.solve(R, A)
    print 'root 1:', A1
    print 'root 2:', A2
    print 'series, root 1:', A1.series(p, 0, 4)
    print 'series, root 2:', A2.series(p, 0, 4)
    # root 1 is always < -1 and will oscillate and grow!

    print
    # Compare A with exact expression exp(-p)
    error_A1 = sm.exp(-p).series(p, 0, 4) - A1.series(p, 0, 4)
    #error_A1 = sm.simplify(error_A1)
    error_A2 = sm.exp(-p).series(p, 0, 4) - A2.series(p, 0, 4)
    #error_A2 = sm.simplify(error_A2)
    print 'Error in A1:', error_A1
    print 'Error in A2:', error_A2

    print
    # Exact Leapfrog solution
    C1, C2 = sm.symbols('C1 C2')
    # Initial condition: u(0)=1 => C1*A**0 + C2*A**0 = 1
    # Need one more condition, say u(n=1) = 1*exp(-p*1)
    q = sm.solve([C1 + C2 - 1, C1*A1 + C2*A2 - 1*sm.exp(-p*1)], [C1, C2])
    print 'Exact Leapfrog solution:'
    solution = q[C1]*A1 + q[C2]*A2
    print solution
    print sm.printing.latex(solution)


def leapfrog_filtered(I, a, b, T, dt, first_step='ForwardEuler', gamma=0.6):
    dt = float(dt)
    Nt = int(round(T/dt))  # no of intervals
    u = np.zeros(Nt+1)
    t = np.linspace(0, Nt*dt, Nt+1)

    n = 0
    u[n] = I
    if first_step == 'ForwardEuler':
        u[n+1] = (1 - dt*a(t[n]))*u[n] + dt*b(t[n])
    elif first_step == 'CrankNicolson':
        u[n+1] = ((1 - 0.5*a(t[n])*dt)*u[n] + \
                  0.5*dt*(b(t[n]) + b(t[n+1])))/\
                  (1 + 0.5*dt*a(t[n+1]))
    else:
        raise ValueError('Wrong first_step=%s' % first_step)

    for n in range(1, Nt):
        u[n+1] = u[n-1] - 2*dt*a(t[n])*u[n] + 2*dt*b(t[n])
        u[n] = u[n] + gamma*(u[n-1] - 2*u[n] + u[n+1])
    return u, t

def demonstrate_cured_instability():
    """Run larger and larger dt and demonstrate small, growing instability."""
    # As demonstrate_instability, but calling leapfrog_filtered

    a = lambda t: 1.0
    b = lambda t: 0
    I = 1.0
    # Solution: exp(-t)
    T = 8
    # Concentrate on coarse meshes
    dt = 0.128
    for k in range(4):
        u, t = leapfrog_filtered(I, a, b, T, dt, first_step='ForwardEuler')
        visualize(u, t, title='timestep: %g' % dt)
        dt *= 2
    # Instability for the last one: A2 root

    # Check longer time integrations
    dt = 0.128
    T = 35
    u, t = leapfrog_filtered(I, a, b, T, dt, first_step='ForwardEuler')
    visualize(u, t, title='timestep: %g' % dt)


if __name__ == '__main__':
    first_test_plot()
    I, c = sm.symbols('I c')
    #analyze(lambda t: c*t + I)
    analyze(lambda t: t**2 + c*t + I)
    #test_linear_solution()
    #convergence_rates()
    #emonstrate_instability()
    #exact_discrete_solution()
    #demonstrate_cured_instability()
