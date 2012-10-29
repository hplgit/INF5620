"""
Approximation of functions by linear combination of basis functions in
function spaces and the least squares method or the collocation method
for determining the coefficients. 2D version.
"""
import sympy as sm
from scitools.std import surfc, savefig, linspace, title, figure, \
     hot, axis, ndgrid, subplot


def least_squares(f, phi, Omega):
    """
    Given a function f(x,y) on a rectangular domain
    Omega=[[xmin,xmax],[ymin,ymax]],
    return the best approximation to f(x,y) in the space V
    spanned by the functions in the list phi.
    """
    N = len(phi) - 1
    A = sm.zeros((N+1, N+1))
    b = sm.zeros((N+1, 1))
    x, y = sm.symbols('x y')
    print '...evaluating matrix...'
    for i in range(N+1):
        for j in range(i, N+1):
            print '(%d,%d)' % (i, j)

            integrand = phi[i]*phi[j]
            I = sm.integrate(integrand,
                             (x, Omega[0][0], Omega[0][1]),
                             (y, Omega[1][0], Omega[1][1]))
            if isinstance(I, sm.Integral):
                # Could not integrate symbolically, use numerical int.
                print 'numerical integration of', integrand
                integrand = sm.lambdify([x,y], integrand)
                I = sm.mpmath.quad(integrand,
                                   [Omega[0][0], Omega[0][1]],
                                   [Omega[1][0], Omega[1][1]])
            A[i,j] = A[j,i] = I
        integrand = phi[i]*f
        I = sm.integrate(integrand,
                         (x, Omega[0][0], Omega[0][1]),
                         (y, Omega[1][0], Omega[1][1]))
        if isinstance(I, sm.Integral):
            # Could not integrate symbolically, use numerical int.
            print 'numerical integration of', integrand
            integrand = sm.lambdify([x,y], integrand)
            I = sm.mpmath.quad(integrand,
                               [Omega[0][0], Omega[0][1]],
                               [Omega[1][0], Omega[1][1]])
        b[i,0] = I
    print
    print 'A:\n', A, '\nb:\n', b
    c = A.LUsolve(b)
    print 'coeff:', c
    # Alternative:
    u = sum(c[i,0]*phi[i] for i in range(len(phi)))
    print 'approximation:', u
    print 'f:', sm.expand(f)
    print sm.latex(A, mode='plain')
    print sm.latex(b, mode='plain')
    print sm.latex([c[i,0] for i in range(len(phi))], mode='plain')
    return u


def comparison_plot(f, u, Omega, plotfile='tmp'):
    """Compare f(x,y) and u(x,y) for x,y in Omega in a plot."""
    x, y = sm.symbols('x y')

    f = sm.lambdify([x,y], f, modules="numpy")
    u = sm.lambdify([x,y], u, modules="numpy")
    # When doing symbolics, Omega can easily contain symbolic expressions,
    # assume .evalf() will work in that case to obtain numerical
    # expressions, which then must be converted to float before calling
    # linspace below
    for r in range(2):
        for s in range(2):
            if not isinstance(Omega[r][s], (int,float)):
                Omega[r][s] = float(Omega[r][s].evalf())

    resolution = 41  # no of points in plot
    xcoor = linspace(Omega[0][0], Omega[0][1], resolution)
    ycoor = linspace(Omega[1][0], Omega[1][1], resolution)
    xv, yv = ndgrid(xcoor, ycoor)
    # Vectorized functions expressions does not work with
    # lambdify'ed functions without the modules="numpy"
    exact  = f(xv, yv)
    approx = u(xv, yv)
    figure()
    surfc(xv, yv, exact, title='f(x,y)',
          colorbar=True, colormap=hot(), shading='flat')
    if plotfile:
        savefig('%s_f.pdf' % plotfile, color=True)
        savefig('%s_f.png' % plotfile)
    figure()
    surfc(xv, yv, approx, title='f(x,y)',
          colorbar=True, colormap=hot(), shading='flat')
    if plotfile:
        savefig('%s_u.pdf' % plotfile, color=True)
        savefig('%s_u.png' % plotfile)

if __name__ == '__main__':
    print 'Module file not meant for execution.'


