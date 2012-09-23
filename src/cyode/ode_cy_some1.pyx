cimport numpy as np
from scitools.std import plot, figure

# Cython does not like lambda function:
# lambda u, t: np.asarray(f_(u, t))
# Made a class wrapper instead

class Wrap:
    def __init__(self, f):
        self.f = f
    def __call__(self, u, t):
        return np.asarray(self.f(u, t))

cpdef RK2(f_, U0, double t):
    f = Wrap(f_)
    cdef np.dnarray t = np.asarray(t)
    cdef int N = len(t)-1
    if isinstance(U0, (list,tuple,np.ndarray)):
        u = np.zeros((N+1, len(U0)))
        u[0] = np.asarray(U0)
    else:
        u = np.zeros(N+1)
        u[0] = U0
             
    cdef int n
    cdef double dt       
    for n in range(N):
        dt = t[n+1] - t[n]
        K1 = dt*f(u[n], t[n])
        K2 = dt*f(u[n] + 0.5*K1, t[n] + 0.5*dt)
        u[n+1] = u[n] + K2
    return u, t

def problem(u, t):
    return [u[1], -u[0]]

class ProblemOpt:
    def __init__(self):
        self.f_return = np.zeros(2)  # avoid remaking numpy arrays

    def __call__(self, u, t):
        self.f_return[0] = u[1]
        self.f_return[1] = -u[0]
        return self.f_return
    
def case(nperiods=4, showplot=False, ftype='class'):
    U0 = [1, 0]
    if ftype == 'class':
        f = ProblemOpt()
    else:
        f = problem
    t0 = time.clock()
    u, t = RK2(f, U0, t=np.linspace(0, nperiods*np.pi, nperiods*30+1))
    t1 = time.clock()
    if showplot:
        plot(t[-200:], u[-200:,0])
        figure()
        plot(u[-200:,0], u[-200:,1])
    return t1-t0

if __name__ == '__main__':
    import time
    cpu = case(1000, showplot=True, ftype='class')
    print 'CPU-time, class version   :', cpu
    cpu = case(1000, showplot=True, ftype='func')
    print 'CPU-time, function version:', cpu



    
