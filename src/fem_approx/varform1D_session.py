import sympy as sm

# Computing with Neumann and Dirichlet conditions: -u''=2
x, C, D = sm.symbols('x C D')
i, j = sm.symbols('i j', integer=True)

integrand = (i+1)*(j+1)*(1-x)**(i+j)
A_ij = sm.integrate(integrand, (x, 0, 1))
A_ij = sm.simplify(A_ij)
print A_ij
psi_i = (1-x)**(i+1)
integrand = 2*psi_i - D*(i+1)*(1-x)**i
b_i = sm.integrate(integrand, (x, 0, 1)) - C*psi_i.subs(x, 0)
b_i = sm.factor(sm.simplify(b_i))
print b_i
print sm.expand(2 - (2+i)*(D+C))

