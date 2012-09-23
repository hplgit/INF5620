from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import sys
problem = sys.argv[1]
version = sys.argv[2]
if not (problem.isdigit() and version.isdigit()):
    print 'Usage: setup.py problem version'
    sys.exit(1)
cymodule = 'ode%s_cy%s' % (problem, version)
# Do not interfere with distutils command-line options
del sys.argv[1]  
del sys.argv[1]

setup(
  name='ODE test',
  ext_modules=[Extension(cymodule, [cymodule + '.pyx'],
                         #libraries=['m'],  # needed for libc.math - no!
                         )],
  cmdclass={'build_ext': build_ext},
)
