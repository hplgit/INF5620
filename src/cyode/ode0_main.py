import sys, time
import numpy as np
import time

version = sys.argv[1]
if int(version) > 12:
    print '1st command-line arg must be version number, not', version
    sys.exit(1)
if version == '0':
    import ode0  as ode # pure Python
else:
    name = 'ode0_cy' + version
    try:
        exec('import ' + name + ' as ode')
        print 'Imported', name
    except ImportError:
        print 'No Cython module', name
        sys.exit(1)
        
def case(N, problem):
    I = 0
    t0 = time.clock()
    T = 5
    dt = float(T)/(N-1)
    u, t = ode.solver(problem, I, dt, T, ode.RK2)
    t1 = time.clock()
    #print 'Final value: u(%g)=%g' % (t, u)
    return t1-t0

try:
    N = int(sys.argv[2])
except IndexError:
    N = 2000001
try:
    problem = sys.argv[3]
except IndexError:
    problem = '1'
cpu = case(N, eval('ode.problem' + problem))
if version == '0':
    module = 'pure Python, ode1.py:'
else:
    module = 'Cython, ode0_cy%s.py:' % version
print 'CPU-time, %s' % module, cpu, 's'
