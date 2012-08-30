# Make a simpler version first

class Problem:
    def __init__():
        self.T, self.I, self.a = 10, 1, 1

class Solver:
    def __init__(self, problem, dt=0.1, theta=0.5):
        self.problem = problem
        self.dt, self.theta = float(dt), theta

    def solve(self):
        N = int(round(problem.T/self.dt))
        self.t = np.linspace(0, problem.T, N+1)

        ...

class Parameters:
    _legal_prms = []

    def __init__(self):
        self.prms = set.default()

    def set(self, **parameters):
        for name in parameters:
            if name in Problem._legal_prms:
                self.prms[name] = parameters[name]

    def get(self, name):
        return self.prms[name]

    def define_command_line_args(self, parser):
        ...

    def load_from_command_line(self, parser):
        for name in Problem._legal_prms:
            self.prms[name] = getattr(parser, name)


class Problem(Parameters):
    """
    Physical parameters for the problem u'=-a*u, u(0)=I,
    with t in [0,T].
    """
    _legal_prms = ['I', 'a', 'T']

    def set_default(self):
        return dict(I=1, a=1, T=10)

class Solver(Parameters):
    _legal_prms = ['dt', 'theta']

    def __init__(self, problem):
        self.problem = problem
        Parameters.__init__()

    def set_default(self):
        return dict(dt=0.1, theta=0.5)

    def solve(self):
        self.t = np.linspace(0, problem.get('T')
