"""
Approximation of functions by linear combination of basis functions in
function spaces and the least squares method or the collocation method
for determining the coefficients.
"""
import sympy as sm
from scitools.std import plot, hold, legend, savefig, linspace, \
     title, xlabel, axis

def least_squares(f, phi, Omega):
    """
    Given a function f(x) on an interval Omega (2-list)
    return the best approximation to f(x) in the space V
    spanned by the functions in the list phi.
    """
    N = len(phi) - 1
    A = sm.zeros((N+1, N+1))
    b = sm.zeros((N+1, 1))
    x = sm.Symbol('x')
    print '...evaluating matrix...',
    for i in range(N+1):
        for j in range(i, N+1):
            print '(%d,%d)' % (i, j)

            integrand = phi[i]*phi[j]
            I = sm.integrate(integrand, (x, Omega[0], Omega[1]))
            if isinstance(I, sm.Integral):
                # Could not integrate symbolically, use numerical int.
                print 'numerical integration of', integrand
                integrand = sm.lambdify([x], integrand)
                I = sm.mpmath.quad(integrand, [Omega[0], Omega[1]])
            A[i,j] = A[j,i] = I
        integrand = phi[i]*f
        I = sm.integrate(integrand, (x, Omega[0], Omega[1]))
        if isinstance(I, sm.Integral):
            # Could not integrate symbolically, use numerical int.
            print 'numerical integration of', integrand
            integrand = sm.lambdify([x], integrand)
            I = sm.mpmath.quad(integrand, [Omega[0], Omega[1]])
        b[i,0] = I
    print
    print 'A:\n', A, '\nb:\n', b
    c = A.LUsolve(b)  # symbolic solve
    print 'coeff:', c

    # c is a sympy Matrix object, numbers are in c[i,0]
    u = sum(c[i,0]*phi[i] for i in range(len(phi)))
    print 'approximation:', u
    return u

def numerical_linsys_solve(A, b, floating_point_calc='sumpy'):
    """
    Given a linear system Au=b as sympy arrays, solve the
    system using different floating-point software.
    floating_point_calc may be 'sympy', 'numpy.float64',
    'numpy.float32'.
    This function is used to investigate ill-conditioning
    of linear systems arising from approximation methods.
    """
    if floating_point_calc == 'sympy':
        #sm.mpmath.mp.dsp = 10  # does not affect the computations here
        A = sm.mpmath.fp.matrix(A)
        b = sm.mpmath.fp.matrix(b)
        print 'A:\n', A, '\nb:\n', b
        c = sm.mpmath.fp.lu_solve(A, b)
        #c = sm.mpmath.lu_solve(A, b) # more accurate
        print 'sympy.mpmath.fp.lu_solve:', c
    elif floating_point_calc.startswith('numpy'):
        import numpy as np
        # Double precision (float64) by default
        A = np.array(A.evalf())
        b = np.array(b.evalf())
        if floating_point_calc == 'numpy.float32':
            # Single precision
            A = A.astype(np.float32)
            b = b.astype(np.float32)
        c = np.linalg.solve(A, b)
        print 'numpy.linalg.solve, %s:' % floating_point_calc, c


def least_squares_orth(f, phi, Omega):
    """
    Same as least_squares, but for orthogonal
    basis such that one avoids calling up standard
    Gaussian elimination.
    """
    N = len(phi) - 1
    A = [0]*(N+1)       # plain list to hold symbolic expressions
    b = [0]*(N+1)
    x = sm.Symbol('x')
    print '...evaluating matrix...',
    for i in range(N+1):
        print '(%d,%d)' % (i, i)
        A[i] = sm.integrate(phi[i]**2, (x, Omega[0], Omega[1]))

        # Fallback on numerical integration if f*phi is too difficult
        # to integrate
        integrand = phi[i]*f
        I = sm.integrate(integrand,  (x, Omega[0], Omega[1]))
        if isinstance(I, sm.Integral):
            print 'numerical integration of', integrand
            integrand = sm.lambdify([x], integrand)
            I = sm.mpmath.quad(integrand, [Omega[0], Omega[1]])
        b[i] = I
    print 'A:\n', A, '\nb:\n', b
    c = [b[i]/A[i] for i in range(len(b))]
    print 'coeff:', c
    u = 0
    for i in range(len(phi)):
        u += c[i]*phi[i]
    # Alternative:
    # u = sum(c[i,0]*phi[i] for i in range(len(phi)))
    print 'approximation:', u
    return u

def interpolation(f, phi, points):
    """
    Given a function f(x), return the approximation to
    f(x) in the space V, spanned by phi, that interpolates
    f at the given points. Must have len(points) = len(phi)
    """
    N = len(phi) - 1
    A = sm.zeros((N+1, N+1))
    b = sm.zeros((N+1, 1))
    # Wrap phi and f in Python functions rather than expressions
    # so that we can evaluate phi at points[i] (alternative to subs?)
    x = sm.Symbol('x')
    phi = [sm.lambdify([x], phi[i]) for i in range(N+1)]
    f = sm.lambdify([x], f)
    print '...evaluating matrix...'
    for i in range(N+1):
        for j in range(N+1):
            print '(%d,%d)' % (i, j)
            A[i,j] = phi[j](points[i])
        b[i,0] = f(points[i])
    print
    print 'A:\n', A, '\nb:\n', b
    c = A.LUsolve(b)
    print 'coeff:', c
    u = 0
    for i in range(len(phi)):
        u += c[i,0]*phi[i](x)
    # Alternative:
    # u = sum(c[i,0]*phi[i] for i in range(len(phi)))
    print 'approximation:', u
    return u

collocation = interpolation  # synonym in this module

def comparison_plot(f, u, Omega, filename='tmp.pdf',
                    plot_title='', ymin=None, ymax=None,
                    u_legend='approximation'):
    """Compare f(x) and u(x) for x in Omega in a plot."""
    x = sm.Symbol('x')
    print 'f:', f

    f = sm.lambdify([x], f, modules="numpy")
    u = sm.lambdify([x], u, modules="numpy")
    if len(Omega) != 2:
        raise ValueError('Omega=%s must be an interval (2-list)' % str(Omega))
    # When doing symbolics, Omega can easily contain symbolic expressions,
    # assume .evalf() will work in that case to obtain numerical
    # expressions, which then must be converted to float before calling
    # linspace below
    if not isinstance(Omega[0], (int,float)):
        Omega[0] = float(Omega[0].evalf())
    if not isinstance(Omega[1], (int,float)):
        Omega[1] = float(Omega[1].evalf())

    resolution = 401  # no of points in plot
    xcoor = linspace(Omega[0], Omega[1], resolution)
    # Vectorized functions expressions does not work with
    # lambdify'ed functions without the modules="numpy"
    exact  = f(xcoor)
    approx = u(xcoor)
    plot(xcoor, approx, '-')
    hold('on')
    plot(xcoor, exact, '-')
    legend([u_legend, 'exact'])
    title(plot_title)
    xlabel('x')
    if ymin is not None and ymax is not None:
        axis([xcoor[0], xcoor[-1], ymin, ymax])
    savefig(filename)

if __name__ == '__main__':
    print 'Module file not meant for execution.'


