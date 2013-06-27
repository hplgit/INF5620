from numpy import *
#from matplotlib.pyplot import *
from scitools.easyviz.matplotlib_ import *

def A_exact(C, p):
    return exp(-C*p**2)

def A_FE(C, p, FEM=False):
    f = 1 + (2./3)*sin(p/2)**2 if FEM else 1.
    return 1 - 4*C*sin(p/2)**2/f

def A_BE(C, p, FEM=False):
    f = 1 + (2./3)*sin(p/2)**2 if FEM else 1.
    return 1/(1 + 4*C*sin(p/2)**2/f)

def A_CN(C, p, FEM=False):
    f = 1 + (2./3)*sin(p/2)**2 if FEM else 1.
    return (1 - 2*C*sin(p/2)**2/f)/(1 + 2*C*sin(p/2)**2/f)

def compare_plot(C, p, FEM=False):
    figure()
    plot(p, A_BE(C, p, FEM),
         p, A_exact(C, p),
         p, A_CN(C, p, FEM),
         p, A_FE(C, p, FEM),)
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
    text = '_FEM' if FEM else '_FDM'
    savefig('A_C%s%s.pdf' % (str(C).replace('.', ''), text))
    savefig('A_C%s%s.png' % (str(C).replace('.', ''), text))

def compare_FEM_FDM(C, p):
    for method in 'FE', 'BE', 'CN':
        figure()
        plot(p, eval('A_%s(C, p, FEM=True)' % method),
             p, eval('A_%s(C, p, FEM=False)' % method),
             p, A_exact(C, p),)
        legend([method + '-FEM', method + '-FDM', 'exact'])
        title('C=%g' % C)
        axis([p[0], p[-1], -1.2, 1])
        xlabel('$p=k\Delta x$')
        savefig('A_C%s_%s_FEMvsFDM.pdf' % (str(C).replace('.', ''), method))
        savefig('A_C%s_%s_FEMvsFDM.png' % (str(C).replace('.', ''), method))

p = linspace(0, pi, 101)
for C in 20, 2, 0.5, 0.25, 0.1, 0.01:
    compare_plot(C, p, FEM=False)
    compare_plot(C, p, FEM=True)
    #compare_FEM_FDM(C, p)


FEM = False
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
doconce combine_images A_C20_FDM.pdf A_C2_FDM.pdf diffusion_A_C20_C2_FDM.pdf
doconce combine_images A_C20_FDM.png A_C2_FDM.png diffusion_A_C20_C2_FDM.png

doconce combine_images A_C05_FDM.png A_C025_FDM.png diffusion_A_C05_C025_FDM.png
doconce combine_images A_C05_FDM.pdf A_C025_FDM.pdf diffusion_A_C05_C025_FDM.pdf

doconce combine_images A_C01_FDM.pdf A_C001_FDM.pdf diffusion_A_C01_C001_FDM.pdf
doconce combine_images A_C01_FDM.png A_C001_FDM.png diffusion_A_C01_C001_FDM.png

doconce combine_images A_C20_FEM.pdf A_C2_FEM.pdf diffusion_A_C20_C2_FEM.pdf
doconce combine_images A_C20_FEM.png A_C2_FEM.png diffusion_A_C20_C2_FEM.png

doconce combine_images A_C05_FEM.png A_C025_FEM.png diffusion_A_C05_C025_FEM.png
doconce combine_images A_C05_FEM.pdf A_C025_FEM.pdf diffusion_A_C05_C025_FEM.pdf

doconce combine_images A_C01_FEM.pdf A_C001_FEM.pdf diffusion_A_C01_C001_FEM.pdf
doconce combine_images A_C01_FEM.png A_C001_FEM.png diffusion_A_C01_C001_FEM.png
"""
