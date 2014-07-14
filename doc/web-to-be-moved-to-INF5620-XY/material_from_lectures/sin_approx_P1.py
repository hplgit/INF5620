import sympy as sm
"""
Exercise 13: Perform symbolic finite element computations
Perform hand calculations to find formulas for the coefficient matrix and right-
hand side when approximating f(x) = sin(x) on = [0;pi] by two P1 elements
of size pi/2. Solve the system and compare u(pi/2) with the exact value 1.
"""

x = sm.Symbol('x')
A = sm.zeros((3,3))
f = sm.sin

phi_0 = 1 - (2*x)/sm.pi
phi_1l = 2*x/sm.pi
phi_1r = 2 - (2*x)/sm.pi
phi_2 = x/(sm.pi/2) - 1
node_0 = 0
node_1 = sm.pi/2
node_2 = sm.pi

A[0,1]=sm.integrate(phi_0*phi_1l, (x, node_0, node_1))
A[1,0]=A[0,1]

A[1,2]= sm.integrate(phi_1r*phi_2, (x, node_1, node_2))
A[2,1]=A[1,2]

A[0,0]= sm.integrate(phi_0**2,  (x, node_0, node_1))
A[1,1]= sm.integrate(phi_1l**2, (x, node_0, node_1)) + \
        sm.integrate(phi_1r**2, (x, node_1, node_2))
A[2,2]= sm.integrate(phi_2**2,  (x, node_1, node_2))

print 'A:\n', A  # Compare with general matrix, h=pi/2

b = sm.zeros((3,1))

b[0] = sm.integrate(phi_0*f(x),  (x, node_0, node_1))
b[1] = sm.integrate(phi_1l*f(x), (x, node_0, node_1)) + \
       sm.integrate(phi_1r*f(x), (x, node_1, node_2))
b[2] = sm.integrate(phi_2*f(x),  (x, node_1, node_2))

print 'b:\n', b

c = A.LUsolve(b)
print 'c:\n', c

for i in range(len(c)):
    print 'c[%d]=%g' % (i, c[i].evalf())
print 'u(pi/2)=c[1]'
