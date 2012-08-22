import odespy
from numpy import *
from matplotlib.pyplot import *

def f(u, t):
    return -a*u

def exact_solution(t):
    return I*exp(-a*t)

I = 1; a = 2
methods = [('ThetaRule', dict(theta=0.5)),
           'RK2', 'RKFehlberg',]
T = 5
dt = 0.5
N = int(round(T/dt))
t = np.linspace(0, T, N+1)
t_fine = np.linspace(0, T, 10001)
legends = []

for method in methods:
    if isinstance(method, (list,tuple)):
        method, parameters = method
    else:
        parameters = {}
    solver = eval('odespy.' + method + '(f, **parameters)')
    solver.set_initial_condition(I)

    def terminate(u, t, step_no):
        """Terminate solver when u is close enough to 0."""
        return abs(u[step_no]) < 1E-6

    u, t = solver.solve(t, terminate)
    legends.append(str(solver))
    if method == 'RKFehlberg':
        plot(solver.t_all, solver.u_all, 'ko')
    else:
        plot(t, u)
    hold('on')
    savefig('tmp_odespy.png')
plot(t_fine, exact_solution(t_fine))
legends.append('exact')
legend(legends)
show()

