from fe_approx1D import approximate

from sympy import Symbol, sin, tanh, pi
x = Symbol('x')

def approx(functionclass='intro'):
    """
    Exemplify approximating various functions by various choices
    of finite element functions
    """
    if functionclass == 'intro':
        approximate(x*(1-x), symbolic=True,
                    d=1, n_e=2, filename='fe_p1_x2_2e')
        approximate(x*(1-x), symbolic=True, numint='Trapezoidal',
                    d=1, n_e=2, filename='fe_p1_x2_2eT')
        approximate(x*(1-x), symbolic=True,
                    d=2, n_e=4, filename='fe_p1_x2_2e')
        approximate(x*(1-x), symbolic=True, numint='Simpson',
                    d=2, n_e=4, filename='fe_p1_x2_2eS')
        #sys.exit(1)
        approximate(x*(1-x), symbolic=False,
                    d=1, n_e=2, filename='fe_p1_x2_2e')
        approximate(x*(1-x), symbolic=False,
                    d=1, n_e=8, filename='fe_p1_x2_8e')

    elif functionclass == 'special':
        # Does not work well because Heaviside cannot be analytically
        # integrated (not important) and numpy cannot evaluate
        # Heaviside (serious for plotting)
        from sympy import Heaviside, Rational
        approximate(Heaviside(x - Rational(1,2)), symbolic=False,
                    d=1, n_e=2, filename='fe_p1_H_2e')
        approximate(Heaviside(x - Rational(1,2)), symbolic=False,
                    d=1, n_e=4, filename='fe_p1_H_4e')
        approximate(Heaviside(x - Rational(1,2)), symbolic=False,
                    d=1, n_e=8, filename='fe_p1_H_8e')
    elif functionclass == 'easy':
        approximate(1-x,        symbolic=False,
                    d=1, n_e=4, filename='fe_p1_x_4e')
        approximate(x*(1-x),    symbolic=False,
                    d=1, n_e=4, filename='fe_p1_x2_4e')
        approximate(x*(1-x),    symbolic=False,
                    d=2, n_e=4, filename='fe_p2_x2_4e')
        approximate(x*(1-x)**8, symbolic=False,
                    d=1, n_e=4, filename='fe_p1_x9_4e')
        approximate(x*(1-x)**8, symbolic=False,
                    d=1, n_e=8, filename='fe_p1_x9_8e')
        approximate(x*(1-x)**8, symbolic=False,
                    d=2, n_e=2, filename='fe_p2_x9_2e')
        approximate(x*(1-x)**8, symbolic=False,
                    d=2, n_e=4, filename='fe_p2_x9_4e')
    elif functionclass == 'hard':
        # Note: it takes time to approximate these cases...
        """
        approximate(sin(pi*x),  symbolic=False,
                    d=1, n_e=4, filename='fe_p1_sin_4e')
        approximate(sin(pi*x),  symbolic=False,
                    d=2, n_e=2, filename='fe_p2_sin_2e')
        approximate(sin(pi*x),  symbolic=False,
                    d=2, n_e=4, filename='fe_p2_sin_4e')
        """
        approximate(sin(pi*x)**2,     symbolic=False,
                    d=1, n_e=4, filename='fe_p1_sin2_4e')
        approximate(sin(pi*x)**2,     symbolic=False,
                    d=2, n_e=2, filename='fe_p2_sin2_2e')
        approximate(sin(pi*x)**2,     symbolic=False,
                    d=2, n_e=4, filename='fe_p2_sin2_4e')
import sys
try:
    task = sys.argv[1]
except IndexError:
    print 'Usage: %s approx_intro | approx_easy | approx_hard' % sys.argv[0]
    sys.exit(1)

approx(task.split('_')[1])


