    A = A.evalf()
    b = b.evalf()
    print 'A:\n', A, '\nb:\n', b
    print 'dps:', sm.mpmath.mp.dps
    #c = A.LUsolve(b)
    c = sm.mpmath.lu_solve(A, b)


def line(p0, p1, x):
    """
    Return y(x) = a*x + b, where a and b are determined such that
    y(x) goes through the points p0 and p1:
    y(p0[0])=p0[1] and y(p0[0])=p0[1].
    """
    x0, y0 = p0
    x1, y1 = p1
    return y0 + ((y1 - y0)/(x1 - x0))*(x - x0)
    
def phi_hat(x, i, nodes):
    """
    Evaluate global linear finite element function, associated
    with node i.
    """
    N = len(nodes)-1
    print i, x
    if i > 0 and nodes[i-1] <= x <= nodes[i]:
        return line((nodes[i-1], 0), (nodes[i], 1), x)
    elif i < N and nodes[i] <= x <= nodes[i+1]:
        return line((nodes[i], 1), (nodes[i+1], 0), x)
    else:
        return 0.


def plot_phi_hat_linear(i, nodes):
    n = 3
    x = concatenate([linspace(nodes[j], nodes[j+1], n) \
                     for j in range(len(nodes)-1)])
    y = [phi_hat(_x, i, nodes) for _x in x]
    plot(x, y)

def plot_some_phis():
    nodes = linspace(0, 2, 9)
    for i in (0, 2, 3, len(nodes)-1):
        plot_phi_hat_linear(i, nodes)
        hold('on')
        legend(int(i))

def phi_symbolic():
    import sympy as sp
    x, x_i, x_ip1, x_im1, h = sp.symbols('x x_i x_ip1 x_im1 h')
    phi_L = line((x_im1,sp.Real(0)), (x_i,sp.Real(1)), x)
    phi_R = line((x_i,sp.Real(1)), (x_ip1,sp.Real(0)), x)
    phi_i = Piecewise((0, x < x_im1), (phi_L, 0 <= x <= x_i),
                      (phi_R, x_i <= x <= x_ip1), (0, x > x_ip1))
    print phi_L
    print phi_R
    print phi_i

#phi_symbolic()

def phi_hat_ref(r, X):
    if r == 0:
        return (1 - X)/2
    elif r == 1:
        return (1 + X)/2
    else:
        raise IndexError('r=%d is not 0 or 1' % r)


# ----------- Functions for arbitrary polynomial degree ------------

