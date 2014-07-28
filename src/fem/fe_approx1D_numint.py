from fe_approx1D import *
import sys

# Extended versions with numerical integration (Midpoint, Trap., Simpson)
# (note that these functions overwrite those imported above!)

def mesh_uniform(n_e, d, Omega=[0,1]):
    """
    Return a 1D finite element mesh on Omega with n_e elements of
    the polynomial degree d. The nodes are uniformly spaced.
    Return vertices (vertices), local vertex to global
    vertex mapping (cells), and local to global degree of freedom
    mapping (dof_map).
    """
    vertices = np.linspace(Omega[0], Omega[1], n_e + 1).tolist()
    dof_map = [[e*d + i for i in range(d+1)] for e in range(n_e)]
    cells = [[e, e+1] for e in range(n_e)]
    return vertices, cells, dof_map

def mesh_uniform_symbolic(n_e, d, Omega=[0,1]):
    """
    As mesh, but use using symbols for the coordinates
    (rational expressions with h as the uniform element length).
    """
    h = sm.Symbol('h')  # element length
    dx = h*sm.Rational(1, d)  # node spacing
    vertices = [Omega[0] + i*dx for i in range(n_e + 1)]
    dof_map = [[e*d + i for i in range(d+1)] for e in range(n_e)]
    cells = [[e, e+1] for e in range(n_e)]
    return vertices, cells, dof_map

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
            if symbolic:
                I = sm.integrate(f*phi[r]*detJ, (X, -1, 1))
            if not symbolic or isinstance(I, sm.Integral):
                # Ensure h is numerical
                h = Omega_e[1] - Omega_e[0]
                detJ = h/2
                #integrand = sm.lambdify([X], f*phi[r]*detJ, modules='sympy')
                integrand = sm.lambdify([X], f*phi[r]*detJ)
                #integrand = integrand.subs(sm.pi, np.pi)
                # integrand may still contain symbols like sm.pi that
                # prevents numerical evaluation...
                try:
                    I = sm.mpmath.quad(integrand, [-1, 1])
                except Exception as e:
                    print 'Could not integrate f*phi[r] numerically:'
                    print e
                    sys.exit(0)
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


def assemble(vertices, cells, dof_map, phi, f,
             symbolic=True, numint=None):
    import sets
    n_n = len(list(set(np.array(dof_map).ravel())))
    n_e = len(cells)
    if symbolic:
        A = sm.zeros((n_n, n_n))
        b = sm.zeros((n_n, 1))    # note: (n_n, 1) matrix
    else:
        A = np.zeros((n_n, n_n))
        b = np.zeros(n_n)
    for e in range(n_e):
        Omega_e = [vertices[cells[e][0]], vertices[cells[e][1]]]
        A_e = element_matrix(phi[e], Omega_e, symbolic, numint)
        b_e = element_vector(f, phi[e], Omega_e, symbolic, numint)
        #print 'element', e
        #print b_e
        for r in range(len(dof_map[e])):
            for s in range(len(dof_map[e])):
                A[dof_map[e][r],dof_map[e][s]] += A_e[r,s]
            b[dof_map[e][r]] += b_e[r]
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
            print 'Numerical rule %s is not supported' % numint
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
            print 'Numerical rule %s is not supported' % numint
            numint = None


    if symbolic:
        vertices, cells, dof_map = mesh_uniform_symbolic(n_e, d, Omega)
    else:
        vertices, cells, dof_map = mesh_uniform(n_e, d, Omega)

    # phi is a list where phi[e] holds the basis in cell no e
    # (this is required by assemble, which can work with
    # meshes with different types of elements).
    # len(dof_map[e]) is the number of nodes in cell e,
    # and the degree of the polynomial is len(dof_map[e])-1
    phi = [basis(len(dof_map[e])-1) for e in range(n_e)]

    print 'phi basis (reference element):\n', phi
    A, b = assemble(vertices, cells, dof_map, phi, f,
                    symbolic=symbolic, numint=numint)

    print 'cells:', cells
    print 'vertices:', vertices
    print 'dof_map:', dof_map
    print 'A:\n', A
    print 'b:\n', b
    #print sm.latex(A, mode='plain')
    #print sm.latex(b, mode='plain')

    if symbolic:
        c = A.LUsolve(b)
    else:
        c = np.linalg.solve(A, b)

    print 'c:\n', c

    if not symbolic:
        print 'Plain interpolation/collocation:'
        x = sm.Symbol('x')
        f = sm.lambdify([x], f, modules='numpy')
        try:
            f_at_vertices = [f(xc) for xc in vertices]
            print f_at_vertices
        except Exception as e:
            print 'could not evaluate f numerically:'
            print e
    # else: nodes are symbolic so f(nodes[i]) only makes sense
    # in the non-symbolic case

    if not symbolic and filename is not None:
        title = 'P%d, n_e=%d' % (d, n_e)
        if numint is None:
            title += ', exact integration'
        else:
            title += ', integration: %s' % numint
        xf = np.linspace(Omega[0], Omega[1], 10001)
        U = np.asarray(c)
        xu, u = u_glob(U, vertices, cells, dof_map, 51)
        from scitools.std import plot
        plot(xu, u, 'r-',
             xf, f(xf), 'b-',
             legend=('u', 'f'),
             title=title)
        savefig(filename + '.pdf')
        savefig(filename + '.png')
    return c

def u_glob(U, vertices, cells, dof_map,
           resolution_per_element=51):
    """
    Compute (x, y) coordinates of a curve y = u(x), where u is a
    finite element function: u(x) = sum_i of U_i*phi_i(x).
    Method: Run through each element and compute curve coordinates
    over the element.
    """
    x_patches = []
    u_patches = []
    for e in range(len(cells)):
        Omega_e = [vertices[cells[e][0]], vertices[cells[e][1]]]
        local_nodes = dof_map[e]
        d = len(local_nodes) - 1
        X = np.linspace(-1, 1, resolution_per_element)
        x = affine_mapping(X, Omega_e)
        x_patches.append(x)
        u_element = 0
        for r in range(len(local_nodes)):
            i = local_nodes[r]  # global node number
            u_element += U[i]*phi_r(r, X, d)
        u_patches.append(u_element)
    x = np.concatenate(x_patches)
    u = np.concatenate(u_patches)
    return x, u

if __name__ == '__main__':
    import sys
    from scitools.misc import function_UI
    cmd = function_UI(
        [phi_r, u_glob, element_matrix, element_vector,
         exemplify_element_matrix_vector, assemble, approximate],
        sys.argv)
    x = sm.Symbol('x')  # needed in eval when expression f contains x
    eval(cmd)
