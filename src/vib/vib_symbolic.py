from sympy import *
x, w, b = symbols('x w b')
eq = x**2 + b*x + w**2
#s = solve(eq == 0, x)
#u = expand(exp(s[1]), complex=True)
#u = re(u)  # im(u)

# not smart enough for this:
"""
t = Symbol('t', real=True)
w_tilde = Symbol('w_tilde', real=True)
dt = Symbol('dt', real=True)

B = exp(I*w_tilde*t)
Q = (B.subs(t, t+dt) - 2*B + B.subs(t, t-dt))
print Q
Q = expand(Q/B, complex=True)
print Q
Q = simplify(Q)
print Q
"""

dt, w = symbols('dt w')
w_tilde = asin(w*dt/2).series(dt, 0, 4)*2/dt
print w_tilde
