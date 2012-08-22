from numpy import linspace, exp
#from matplotlib.pyplot import *
from scitools.std import *

def A_exact(p):
    return exp(-p)

def A(p, theta):
    return (1-(1-theta)*p)/(1+theta*p)

def A_leapfrog(p):
    LF = - p + sqrt(p**2 + 1)
    # unstable leapfrog mode:
    #LF2 = - p - sqrt(p**2 + 1)
    return LF

def amplification_factor(names):
    curves = {}
    p = linspace(0, 3, 101)
    curves['exact'] = A_exact(p)
    name2theta = dict(FE=0, BE=1, CN=0.5)
    for name in names:
        if name in name2theta:
            curves[name] = A(p, name2theta[name])
        elif name == 'LF':
            curves[name] = A_leapfrog(p)
    for name in names:
        plot(p, curves[name])
        hold('on')
    plot([p[0], p[-1]], [0, 0], '--')  # A=0 line
    title('Amplification factors')
    grid('on')
    legend(names + ['exact'], loc='lower left', fancybox=True)
    xlabel('a*dt')
    ylabel('A')
    savefig('A_factors.png')
    savefig('A_factors.eps')
    show()

if __name__ == '__main__':
    if len(sys.argv) > 1:
        names = sys.argv[1:]
    else: # default
        names = 'FE BE CN'.split()
    amplification_factor(names)
