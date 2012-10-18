# generate interactive demo session running the script
# below as input to scitools file2interactive

from fe_approx1D import *

import matplotlib.pyplot as plt
# u example with two sine basis functions
x = np.linspace(0, 4, 1001)
phi0 = np.sin(2*np.pi/4*x)
phi1 = np.sin(2*np.pi*x)
phi2 = np.sin(2*np.pi*4*x)
u = 4*phi0 - 0.5*phi1 - 0*phi2
plt.plot(x, phi0, 'r-', label=r"$\varphi_0$")
plt.plot(x, phi1, 'g-', label=r"$\varphi_1$")
#plt.plot(x, phi2, label=r"$\varphi_2$")
plt.plot(x, u, 'b-', label=r"u")
plt.legend()
plt.savefig('tmp.pdf')
plt.savefig('tmp.png')

# u example with hat basis functions
plt.figure()
x = [0, 1.5, 2.5, 3.5, 4]
phi = [np.zeros(len(x)) for i in range(len(x)-2)]
for i in range(len(phi)):
    phi[i][i+1] = 1
#u = 5*x*np.exp(-0.25*x**2)*(4-x)
u = [0, 8, 5, 4, 0]
for i in range(len(phi)):
    plt.plot(x, phi[i], 'r-')  #, label=r'$\varphi_%d$' % i)
    plt.text(x[i+1], 1.2, r'$\varphi_%d$' % i)
plt.plot(x, u, 'b-', label='$u$')
plt.legend(loc='upper left')
plt.axis([0, x[-1], 0, 9])
plt.savefig('tmp2.png')
plt.savefig('tmp2.pdf')
# Mark elements
for xi in x[1:-1]:
    plt.plot([xi, xi], [0, 9], 'm--')
# Mark nodes
#plt.plot(x, np.zeros(len(x)), 'ro', markersize=4)
plt.savefig('tmp3.png')
plt.savefig('tmp3.pdf')
plt.show()
import sys
sys.exit(0)


# "Hand"-integration of element matrix and vector
x, x_m, h, X = sm.symbols('x x_m h X')
A_00 = sm.integrate(h/8*(1-X)**2, (X, -1, 1))
print A_00
print sm.latex(A_00, mode='plain')
A_10 = sm.integrate(h/8*(1+X)*(1-X), (X, -1, 1))
print A_10
print sm.latex(A_10, mode='plain')
A_11 = sm.integrate(h/8*(1+X)**2, (X, -1, 1))
print A_11
print sm.latex(A_11, mode='plain')
x = x_m + h/2*X
b_0 = sm.integrate(h/4*x*(1-x)*(1-X), (X, -1, 1))
b_1 = sm.integrate(h/4*x*(1-x)*(1+X), (X, -1, 1))
print b_0
print b_1
print sm.latex(b_0, mode='plain')
print sm.latex(b_1, mode='plain')

phi = basis(d=1)
phi
element_matrix(phi, Omega_e=[0.1, 0.2], symbolic=True)
element_matrix(phi, Omega_e=[0.1, 0.2], symbolic=False)

h, x = sm.symbols('h x')
nodes = [0, h, 2*h]
elements = [[0, 1], [1, 2]]
phi = basis(d=1)
f = x*(1-x)
A, b = assemble(nodes, elements, phi, f, symbolic=True)
A
b
c = A.LUsolve(b)
c
fn = sm.lambdify([x], f)
[fn(xc) for xc in nodes]

# The corresponding numerical computations, as done by sympy and
# still based on symbolic integration, goes as follows:

nodes = [0, 0.5, 1]
elements = [[0, 1], [1, 2]]
phi = basis(d=1)
x = sm.Symbol('x')
f = x*(1-x)
A, b = assemble(nodes, elements, phi, f, symbolic=False)
A
b
c = A.LUsolve(b)
c

d=1; n_e=8; Omega=[0,1]  # 8 linear elements on [0,1]
phi = basis(d)
f = x*(1-x)
nodes, elements = mesh_symbolic(n_e, d, Omega)
A, b = assemble(nodes, elements, phi, f, symbolic=True)
A

