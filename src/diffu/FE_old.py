from scitools.std import *

def u(t, x, alpha, k, dx):
    # alpha is Elements book notation, not diffusion coeff!
    """Exact discrete solution of the diffusion eq. (Forward Euler)."""
    dt = dx**2*alpha
    n = int(t/dt)
    return (1 - 4*alpha*sin((k*dx/2))**2)**n*sin(k*x)

dx = 1./9
C = alpha = 0.4765
k = 3*pi
x = linspace(0,1,41)
t = 0.1
y1 = u(t, x, alpha, k, dx)
plot(x, y1)
print 'Amp error:', C**2*dx**4*k**4/2
hold('on')
dx = 1./59
alpha = 0.4931
y2 = u(t, x, alpha, k, dx)
plot(x, y2)
print 'Amp error:', C**2*dx**4*k**4/2

# Remarkable higher damping with larger dx!

from sympy import *
k, dx, C = symbols('k dx C')
f = 1 - 4*C*sin((k*dx/2))**2
fs = f.series(dx, 0, 4)
fex = exp(-C*k**2*dx**2).series(dx, 0, 4)
print fs
print fex - fs
fs = 1 - C*dx**2*k**2 + C*dx**4*k**4/12
fs = fs**(1/(C*dx**2))  # this gets too complicated
print fs.series(dx, 3, 0)
# Alternative: plot
t = 0.1
K = f**(t/(C*dx**2))
print '---', K
K = K.subs('k', 3*pi)
print K
K = K.subs('C', 0.5)
print K
f = lambdify([dx], K, modules="numpy")
K = exp(-k**2*t).subs('k', 3*pi).subs('C', 0.5)
print K
fex = lambdify([dx], K, modules="numpy")
dx = linspace(0.001, 0.1, 11)
figure()
y1 = f(dx)
y2 = fex(dx)
print y1
print y2
plot(dx, y1, [dx[0], dx[-1]], [y2, y2])


