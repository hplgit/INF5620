from fe_approx1D import *

# Extended versions with numerical integration (Midpoint, Trap., Simpson)
# (note that these functions overwrite those imported above!)

def mesh2(n_e, d, Omega=[0,1]):
    """
    Return a 1D finite element mesh on Omega with n_e elements of
    the polynomial degree d. The nodes are uniformly spaced.
    Return vertices (vertices), local vertex to global
    vertex mapping (cells), and local to global degree of freedom
    mapping (dof_map).
    """
    vertices = np.linspace(Omega[0], Omega[1], n_e + 1).tolist()
    doc_map = [[e*d + i for i in range(d+1)] for e in range(n_e)]
    cells = [[e, e+1] for e in range(n_e)]
    return vertices, cells, dof_map
    # Not yet used

def element_matrix(phi, Omega_e, symbolic=True, numint=None):
    n = len(phi)
    A_e = sm.zeros((n, n))
    X = sm.Symbol('X')
    if symbolic:
        h = sm.Symbol('h')
    else:
        h = Omega_e[1] - Omega_e[0]
    detJ = h/2  # dx/dX
    if numint is None:
        for r in range(n):
            for s in range(r, n):
                A_e[r,s] = sm.integrate(phi[r]*phi[s]*detJ, (X, -1, 1))
                A_e[s,r] = A_e[r,s]
    else:
        #phi = [sm.lambdify([X], phi[r]) for r in range(n)]
        # Do instead a phi_rj = phi[r].subs(X, Xj) to avoid real numbers
        for r in range(n):
            for s in range(r, n):
                for j in range(len(numint[0])):
                    Xj, wj = numint[0][j], numint[1][j]
                    phi_rj = phi[r].subs(X, Xj)
                    phi_sj = phi[s].subs(X, Xj)
                    A_e[r,s] += phi_rj*phi_sj*detJ*wj
                A_e[s,r] = A_e[r,s]
    return A_e

def element_vector(f, phi, Omega_e, symbolic=True, numint=None):
    n = len(phi)
    b_e = sm.zeros((n, 1))
    # Make f a function of X (via f.subs to avoid real numbers from lambdify)
    X = sm.Symbol('X')
    if symbolic:
        h = sm.Symbol('h')
    else:
        h = Omega_e[1] - Omega_e[0]
    x = (Omega_e[0] + Omega_e[1])/2 + h/2*X  # mapping
    f = f.subs('x', x)
    detJ = h/2
    if numint is None:
        for r in range(n):
            I = sm.integrate(f*phi[r]*detJ, (X, -1, 1))
            if isinstance(I, sm.Integral):
                print 'numerical integration of', f*phi[r]*detJ
                # Ensure h is numerical
                h = Omega_e[1] - Omega_e[0]
                detJ = h/2
                integrand = sm.lambdify([X], f*phi[r]*detJ)
                I = sm.mpmath.quad(integrand, [-1, 1])
            b_e[r] = I
    else:
        #phi = [sm.lambdify([X], phi[r]) for r in range(n)]
        # f contains h from the mapping, substitute X with Xj
        # instead of f = sm.lambdify([X], f)
        for r in range(n):
            for j in range(len(numint[0])):
                Xj, wj = numint[0][j], numint[1][j]
                fj = f.subs(X, Xj)
                phi_rj = phi[r].subs(X, Xj)
                b_e[r] += fj*phi_rj*detJ*wj
    return b_e

def exemplify_element_matrix_vector(f, d, symbolic=True, numint=False):
    Omega_e = [0.1, 0.2]
    A_e = element_matrix(phi, Omega_e=Omega_e,
                         symbolic=symbolic, numint=numint)
    integration_msg = """
    Symbolic integration failed, and then numerical integration
    encountered an undefined symbol (because of the symbolic expressions):
    %s"""
    if symbolic:
        h = sm.Symbol('h')
        Omega_e=[1*h, 2*h]
    try:
        b_e = element_vector(f, phi, Omega_e=Omega_e,
                             symbolic=symbolic, numint=numint)
    except NameError as e:
        raise NameError(integration_msg % e)
    print 'Element matrix:\n', A_e
    print 'Element vector:\n', b_e


def assemble(nodes, elements, phi, f, symbolic=True, numint=None):
    n_n, n_e = len(nodes), len(elements)
    zeros = sm.zeros if symbolic else np.zeros
    A = zeros((n_n, n_n))
    b = zeros((n_n, 1))
    for e in range(n_e):
        Omega_e = [nodes[elements[e][0]], nodes[elements[e][-1]]]
        A_e = element_matrix(phi, Omega_e, symbolic, numint)
        b_e = element_vector(f, phi, Omega_e, symbolic, numint)
        #print 'element', e
        #print b_e
        for r in range(len(elements[e])):
            for s in range(len(elements[e])):
                A[elements[e][r],elements[e][s]] += A_e[r,s]
            b[elements[e][r]] += b_e[r]
    return A, b

def approximate(f, symbolic=False, d=1, n_e=4, numint=None,
                Omega=[0, 1], filename='tmp.pdf'):
    if symbolic:
        if numint == 'Trapezoidal':
            numint = [[sm.S(-1), sm.S(1)], [sm.S(1), sm.S(1)]]  # sympy integers
        elif numint == 'Simpson':
            numint = [[sm.S(-1), sm.S(0), sm.S(1)],
                      [sm.Rational(1,3), sm.Rational(4,3), sm.Rational(1,3)]]
        elif numint == 'Midpoint':
            numint = [[sm.S(0)],  [sm.S(2)]]
        elif numint == 'GaussLegendre2':
            numint = [[-1/sm.sqrt(3), 1/sm.sqrt(3)], [sm.S(1), sm.S(1)]]
        elif numint == 'GaussLegendre3':
            numint = [[-sm.sqrt(sm.Rational(3,5)), 0,
                       sm.sqrt(sm.Rational(3,5))],
                      [sm.Rational(5,9), sm.Rational(8,9),
                       sm.Rational(5,9)]]
        else:
            numint = None
    else:
        if numint == 'Trapezoidal':
            numint = [[-1, 1], [1, 1]]
        elif numint == 'Simpson':
            numint = [[-1, 0, 1], [1./3, 4./3, 1./3]]
        elif numint == 'Midpoint':
            numint = [[0], [2]]
        elif numint == 'GaussLegendre2':
            numint = [[-1/sqrt(3), 1/sqrt(3)], [1, 1]]
        elif numint == 'GaussLegendre3':
            numint = [[-sqrt(3./5), 0, sqrt(3./5)],
                      [5./9, 8./9, 5./9]]
        else:
            numint = None

    phi = basis(d)
    print 'phi basis (reference element):\n', phi
    integration_msg = """
    Symbolic integration failed, and then numerical integration
    encountered an undefined symbol (because of the symbolic expressions):
    %s"""

    if symbolic:
        try:
            nodes, elements = mesh_symbolic(n_e, d, Omega)
        except NameError as e:
            raise NameError(integration_msg % e)
    else:
        nodes, elements = mesh(n_e, d, Omega)

    A, b = assemble(nodes, elements, phi, f,
                    symbolic=symbolic, numint=numint)

    print 'nodes:', nodes
    print 'elements:', elements
    print 'A:\n', A
    print 'b:\n', b
    #print sm.latex(A, mode='plain')
    #print sm.latex(b, mode='plain')

    if symbolic:
        c = A.LUsolve(b)
    else:
        c = np.linalg.solve(A, b)

    print 'c:\n', c

    print 'Plain interpolation/collocation:'
    x = sm.Symbol('x')
    f = sm.lambdify([x], f, modules='numpy')
    try:
        f_at_nodes = [f(xc) for xc in nodes]
    except NameError as e:
        raise NameError('numpy does not support special function:\n%s' % e)
    print f_at_nodes
    if not symbolic and filename is not None:
        xf = np.linspace(Omega[0], Omega[1], 10001)
        U = np.asarray(c)
        xu, u = u_glob(U, elements, nodes)
        from scitools.std import plot
        plot(xu, u, 'r-',
             xf, f(xf), 'b-',
             legend=('u', 'f'),
             savefig=filename)


if __name__ == '__main__':
    import sys
    from scitools.misc import function_UI
    cmd = function_UI(
        [phi_r, u_glob, element_matrix, element_vector,
         exemplify_element_matrix_vector, assemble, approximate],
        sys.argv)
    x = sm.Symbol('x')  # needed in eval when expression f contains x
    eval(cmd)
