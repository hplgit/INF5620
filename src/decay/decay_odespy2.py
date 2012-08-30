import odespy
import numpy as np
import decay_theta
#from matplotlib.pyplot import *
from scitools.std import *

def f(u, t):
    return -a*u

def exact_solution(t):
    return I*exp(-a*t)

I = 1; a = 2; T = 4; dt = 0.8
methods = [('BackwardEuler', {'eps_iter': 0.01,
                              'nonlinear_solver': 'Newton',
                              'jac': lambda u,t: -a}),
           ('CrankNicolson', {'eps_iter': 0.01}),
           ('RK4', {}),
           ('DormandPrince', {'atol': 1E-5, 'rtol': 1E-4}),
           ]


N = int(round(T/dt))
t = np.linspace(0, T, N+1)
t_fine = np.linspace(0, T, 10001)
legends = []
errors = []

for method, parameters in methods:
    solver = eval('odespy.' + method + '(f, **parameters)')
    solver.set_initial_condition(I)

    def terminate(u, t, step_no):
        """Terminate solver when u is close enough to 0."""
        return abs(u[step_no]) < 1E-6

    try:
        print 'Running', str(solver)
        u, t = solver.solve(t, terminate)
    except Exception, e:
        print e
        continue

    error = sqrt(dt*np.sum((exact_solution(t) - u)**2))
    errors.append((method, error))

    legends.append(str(solver))
    if hasattr(solver, 'u_all'):
        plot(solver.t_all, solver.u_all, 'ko')
    else:
        plot(t, u)
    hold('on')
    savefig('tmp_odespy.png')
plot(t_fine, exact_solution(t_fine))
legends.append('exact')
legend(legends)
for method, error in errors:
    print '%-20s %.3E' % (method, error)
show()

