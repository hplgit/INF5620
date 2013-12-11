from numpy import *
#from matplotlib.pyplot import *
from scitools.easyviz.matplotlib_ import *

def A_exact(C, p):
    return exp(-4*C*p**2)

def A_FE(C, p):
    return 1 - 4*C*sin(p)**2

def A_BE(C, p):
    return 1/(1 + 4*C*sin(p)**2)

def A_CN(C, p):
    return (1 - 2*C*sin(p)**2)/(1 + 2*C*sin(p)**2)

def compare_plot(C, p):
    figure()
    plot(p, A_BE(C, p),
         p, A_exact(C, p),
         p, A_CN(C, p),
         p, A_FE(C, p),)
    legend(['BE', 'exact', 'CN', 'FE'])
    title('C=%g' % C)
    print 'C:', C
    if 0.2 >= C > 0.02:
        axis([p[0], p[-1], 0.3, 1])
    elif C <= 0.02:
        axis([p[0], p[-1], 0.75, 1])
    else:
        axis([p[0], p[-1], -1.2, 1])
    xlabel('$p=k\Delta x$')
    savefig('A_C%s.pdf' % (str(C).replace('.', '')))
    savefig('A_C%s.png' % (str(C).replace('.', '')))


p = linspace(0, pi/2, 101)
#for C in 20, 2, 0.5, 0.25, 0.1, 0.01:
#    compare_plot(C, p)

from sympy import *
C, p, dx, dt = symbols('C p dx dt')
#A_err_FE = A_FE(C, p)/A_exact(C, p)
A_err_FE = A_exact(C, p) - A_FE(C, p)
A_err_FE = A_FE(C, p)/A_exact(C, p)
#print 'Error in A, FE:', A_err_FE.series(C, 0, 6)
A_err_FE = A_err_FE.subs('C', 'dt/dx**2').subs('p', 'dx')
print 'Error in A, FE:', A_err_FE.series(dt, 0, 6)
print latex(A_err_FE.series(C, 0, 6))
A_err_BE = A_exact(C, p) - A_BE(C, p)
A_err_BE = A_BE(C, p)/A_exact(C, p)
print 'Error in A, BE:', A_err_BE.series(C, 0, 6)
print latex(A_err_BE.series(C, 0, 6))
A_err_CN = A_exact(C, p) - A_CN(C, p)
A_err_CN = A_CN(C, p)/A_exact(C, p)
print 'Error in A, CN:', A_err_CN.series(C, 0, 6)
print latex(A_err_CN.series(C, 0, 6))

raw_input()

show()

"""
doconce combine_images A_C20.pdf A_C2.pdf diffusion_A_C20_C2.pdf
doconce combine_images A_C20.png A_C2.png diffusion_A_C20_C2.png

doconce combine_images A_C05.png A_C025.png diffusion_A_C05_C025.png
doconce combine_images A_C05.pdf A_C025.pdf diffusion_A_C05_C025.pdf

doconce combine_images A_C01.pdf A_C001.pdf diffusion_A_C01_C001.pdf
doconce combine_images A_C01.png A_C001.png diffusion_A_C01_C001.png
"""
