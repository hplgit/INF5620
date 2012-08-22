# This file is aimed at being run as an interactive session via
# scitools file2interactive decay_analysis_symbolic.py

import sys
from sympy import *

p = Symbol('p')
A_e = exp(-p)

# Demo on Taylor polynomials
A_e.series(p, 0, 6)

"""
# NOTE: rsolve can solve recurrence relations:
a, dt, I, n = symbols('a dt I n')
u = Function('u')
f = u(n+1) - u(n) + dt*a*u(n+1)
rsolve(f, u(n), {u(0): I})
# However, 0 is the answer!
# Experimentation shows that we cannot have symbols dt, a in the
# recurrence equation, just n or numbers.
# Even if we worked with scaled equations, dt is in there,
# rsolve cannot be applied.
"""

# Numerical amplification factor
theta = Symbol('theta')
A = (1-(1-theta)*p)/(1+theta*p)

# Interactive session for demonstrating subs
A.subs(theta, 1)  # Backward Euler factor
A.subs(theta, Rational(1,2))  # Crank-Nicolson factor
A.subs(theta, 0).series(p, 0, 4)
A.subs(theta, 1).series(p, 0, 4)
A.subs(theta, Rational(1,2)).series(p, 0, 4)
A_e.series(p, 0, 4)


# Error in amplification factors
half = Rational(1,2)
FE = A_e.series(p, 0, 4) - A.subs(theta, 0).series(p, 0, 4)
BE = A_e.series(p, 0, 4) - A.subs(theta, 1).series(p, 0, 4)
CN = A_e.series(p, 0, 4) - A.subs(theta, half).series(p, 0, 4)
FE
BE
CN

# Ratio of amplification factors
FE = 1 - (A.subs(theta, 0)/A_e).series(p, 0, 4)
BE = 1 - (A.subs(theta, 1)/A_e).series(p, 0, 4)
CN = 1 - (A.subs(theta, half)/A_e).series(p, 0, 4)
FE
BE
CN

print "Error in solution:"
n, a, dt, t, T = symbols('n a dt t T')
u_e = exp(-p*n)
u_n = A**n
FE = u_e.series(p, 0, 4) - u_n.subs(theta, 0).series(p, 0, 4)
print FE
FE = FE.subs('n', 't/dt').subs(p, 'a*dt')
FE = FE.extract_leading_order(dt)[0][0]
FEi = integrate(FE**2, (t, 0, T))
print FEi, type(FEi)
FEi = sqrt(FEi)
print FEi.series(dt, 0, 3)
FEi = simplify(FEi.series(dt, 0, 3).extract_leading_order(dt)[0][0])
print FEi

sys.exit(0)
#BE = u_e.series(p, 4) - u_n.subs(theta, 1).series(p, 4)
#CN = u_e.series(p, 4) - u_n.subs(theta, half).series(p, 4)
FE
BE
CN
simplify(FE)
simplify(BE)
simplify(CN)



def L2_error(local_error, nterms=4):
    """Sum over n (n -> 1/2*N, n^2 -> 1/3*N^2, N = T/dt)."""
    # DOES NOT WORK
    # Replace n by k
    e2 = local_error
    print '1.', e2
    e2 = polynomials.Polynomial(e2)
    e2 = e2.leading_term(e2)
    k = Symbol('k')
    e2 = e2.subs(n, k)
    print '2.', e2
    a = Symbol('a');  dt = Symbol('dt'); T = Symbol('T'); N = Symbol('N')
    e2 = e2**2
    print '3.', e2
    e2 = e2.series(p, 2*nterms)
    print '5.', e2
    e2 = e2.subs(p, dt*a)
    print '6.', e2
    # Perform sum using formulas for sum of k**s
    for s in range(nterms, 0, -1):
        e2 = e2.subs(k**s, Rational(1,s+1)*N**(s+1))
        print '7-%d' % s, e2
    e2 = e2.subs(N, T/dt)
    print '8.', e2
    e2 = e2.series(dt, nterms+1)
    print '9', e2
    e2 = sqrt(e2)
    print '10', e2
    e2 = e2.series(dt, nterms)
    print '11', e2
    return e2


# Leapfrog scheme, amplification factor (exact numerical solution)
A = Symbol('A')
A1, A2 = solve(A**2 + 2*p*A - 1, A)
A1, A2
A1 = simplify(A1); A2 = simplify(A2)
A1, A2
A1.series(p, 6)
A2.series(p, 6)
#print printing.latex(r1)
#print printing.latex(r2)

# leapfrog solution:
C1 = Symbol('C1')
C2 = Symbol('C2')
q = solve([C1 + C2 - 1, C1*A1 + C2*A2 - exp(-p)], [C1, C2])
#print 'q:', q  # big expression
# series approx
q[C1].series(p, 4)  # C1
q[C2].series(p, 4)  # C2
q = solve([C1 + C2 - 1, C1*A1 + C2*A2 - (1-p)], [C1, C2])  # FE
# series approx (FE)
q[C1].series(p, 4) # C1
q[C2].series(p, 4) # C2

