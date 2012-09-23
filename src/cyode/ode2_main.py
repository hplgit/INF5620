import sys, time
import numpy as np
import time

version = sys.argv[1]
if version not in '0123456789':
    print '1st command-line arg must be version number, not', version
    sys.exit(1)
if version == '0':
    import ode2 as ode  # pure Python
else:
    name = 'ode2_cy' + version
    try:
        exec('import ' + name + ' as ode')
        print 'Imported', name
    except ImportError:
        print 'No Cython module', name
        sys.exit(1)
        
def case(nperiods=4, showplot=False):
    U0 = [1, 0]   # Does not work well with float arrays and Cython
    U0 = [1., 0.]
    f = ode.Problem1()
    t0 = time.clock()
    time_points = np.linspace(0, nperiods*2*np.pi, nperiods*30+1)
    u, t = ode.solver(f, U0, time_points, ode.RK2)
    t1 = time.clock()
    if showplot and nperiods < 8:
        from scitools.std import plot, figure
        plot(t[-200:], u[-200:,0])
        figure()
        plot(u[-200:,0], u[-200:,1])
    return t1-t0, time_points.size

try:
    nperiods = int(sys.argv[2])
except IndexError:
    nperiods = 30000
    #nperiods = 1000
cpu, n = case(nperiods, showplot=True)
if version == '0':
    module = 'pure Python, ode.py:'
else:
    module = 'Cython, ode_cy%s.py:' % version
print 'CPU-time, %s' % module, cpu, 's', n
