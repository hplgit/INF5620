import sympy as sm
x, L, C, D, c_0, c_1, = sm.symbols('x L C D c_0 c_1')

f = sm.sin(x)

def model1(f, L, D):
    """Solve -u'' = f(x), u(0)=0, u(L)=D."""
    u_x = - sm.integrate(f, (x, 0, x)) + c_0
    u = sm.integrate(u_x, (x, 0, x)) + c_1
    r = sm.solve([u.subs(x, 0)-0, u.subs(x,L)-D], [c_0, c_1])
    u = u.subs(c_0, r[c_0]).subs(c_1, r[c_1])
    u = sm.simplify(sm.expand(u))
    return u

def model2(f, L, C, D):
    """Solve -u'' = f(x), u'(0)=C, u(L)=D."""
    u_x = - sm.integrate(f, (x, 0, x)) + c_0
    u = sm.integrate(u_x, (x, 0, x)) + c_1
    r = sm.solve([sm.diff(u,x).subs(x, 0)-C, u.subs(x,L)-D], [c_0, c_1])
    u = u.subs(c_0, r[c_0]).subs(c_1, r[c_1])
    u = sm.simplify(sm.expand(u))
    return u

def model3(f, a, L, C, D):
    """Solve -(a*u')' = f(x), u(0)=C, u(L)=D."""
    au_x = - sm.integrate(f, (x, 0, x)) + c_0
    u = sm.integrate(au_x/a, (x, 0, x)) + c_1
    r = sm.solve([u.subs(x, 0)-C, u.subs(x,L)-D], [c_0, c_1])
    u = u.subs(c_0, r[c_0]).subs(c_1, r[c_1])
    u = sm.simplify(sm.expand(u))
    return u


def test1():
    f = 2
    u = model1(f, L, D)
    print 'model1:', u, u.subs(x, 0), u.subs(x, L)
    print sm.latex(u, mode='plain')
    u = model2(f, L, C, D)
    print 'model2:', u, sm.diff(u, x).subs(x, 0), u.subs(x, L)
    print sm.latex(u, mode='plain')
    u = model3(0, 1+x**2, L, C, D)
    print 'model3:', u, u.subs(x, 0), u.subs(x, L)
    print sm.latex(u, mode='plain')

def test2():
    pass

if __name__ == '__main__':
    test1()


