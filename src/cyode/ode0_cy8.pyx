"""
ODE integration restricted to scalar ODEs.
No use of arrays.
Cython version with declaration of variables.
Function objects transferred as arguments are 
made as class instances,
cf. http://docs.cython.org/src/tutorial/cdef_classes.html.
Dropped except -10001 in functions. (20% effect)

Testing Problem2 with exp(t) function and various types of
libraries for such functions. Here: math.h.
"""
cdef class Problem:
    cpdef double rhs(self, double u, double t):
        return 0

cdef class Problem1(Problem):
    cpdef double rhs(self, double u, double t):
        return -u + 1

cdef extern from "math.h":
    double exp(double)

cdef class Problem2(Problem):
    cpdef double rhs(self, double u, double t):
        return - u + exp(-2*t)

cdef class ODEMethod:
    cpdef double advance(self, double u, double t, Problem p,
                        double dt):
        return 0

cdef class Method_RK2(ODEMethod):
    cpdef double advance(self, double u, double t, Problem p,
                         double dt):
        cdef double K1, K2, unew
        K1 = dt*p.rhs(u, t)
        K2 = dt*p.rhs(u + 0.5*K1, t + 0.5*dt)
        unew = u + K2
        return unew
    

cpdef solver(Problem f, double I, double dt, 
             double T, ODEMethod method):
    cdef int N = int(round(float(T)/dt))
    cdef double u = I  # previous time step
    cdef double t = 0
    cdef int n
    for n in xrange(N):
        u = method.advance(u, t, f, dt)
        t += dt
    return u, t


# Create names compatible with ode0.py
RK2 = Method_RK2()
problem1 = Problem1()
problem2 = Problem2()




    
