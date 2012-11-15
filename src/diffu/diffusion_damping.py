from matplotlib.pyplot import *
from numpy import *

def component(Q, x, t, k, a=1):
    return Q*exp(-a*k**2*t)*sin(k*x)

def u(x, t):
    return component(1, x, t, pi, 1) + component(0.1, x, t, 100*pi, 1)

x = linspace(0, 1, 2001)
a = 1
amplitudes = array([0.1, 0.01])
k = 100*pi
times1 = log(amplitudes)/(-a*k**2)
k = pi
times2 = log(amplitudes)/(-a*k**2)
times = [0] + times1.tolist() + times2.tolist()

for t in times:
    figure()
    plot(x, u(x,t))
    title('t=%.2E' % t)
    axis([0, 1, -0.1, 1.1])
    savefig('tmp_%.2E.pdf' % t)
    savefig('tmp_%.2E.png' % t)
show()

# doconce combine_images tmp_0.00E+00.pdf tmp_4.67E-05.pdf tmp_2.33E-01.pdf tmp_4.67E-01.pdf diffusion_damping.pdf
# doconce combine_images tmp_0.00E+00.png tmp_4.67E-05.png tmp_2.33E-01.png tmp_4.67E-01.png diffusion_damping.png



