import sys, time
import numpy as np
import time

version = sys.argv[1]
if version not in '0123456789':
    print '1st command-line arg must be version number, not', version
    sys.exit(1)
if version == '0':
    import ode1 as ode  # pure Python
else:
    name = 'ode1_cy' + version
    try:
        exec('import ' + name + ' as ode')
        print 'Imported', name
    except ImportError:
        print 'No Cython module', name
        sys.exit(1)
        
def case(N, problem=ode.problem2, showplot=False):
    I = 0
    t0 = time.clock()
    t = np.linspace(0, 5, N)
    u, t = ode.solver(problem, I, t, ode.RK2)
    t1 = time.clock()
    if showplot:
        from scitools.std import plot
        plot(t, u)
    return t1-t0

try:
    N = int(sys.argv[2])
except IndexError:
    N = 1000001
cpu = case(N, ode.problem1, showplot=False)
if version == '0':
    module = 'pure Python, ode.py:'
else:
    module = 'Cython, ode1_cy%s.py:' % version
print 'CPU-time, %s' % module, cpu, 's'
