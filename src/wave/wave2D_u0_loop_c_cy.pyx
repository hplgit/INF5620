import numpy as np
cimport numpy as np
cimport cython

cdef extern from "wave2D_u0_loop_c.h":
    void advance(double* u, double* u_1, double* u_2, double* f,
                 double* x, double* y, double* t,
                 double Cx2, double Cy2, double dt2,
                 int Nx, int Ny, int N)

@cython.boundscheck(False)
@cython.wraparound(False)
def advance_cwrap(
    np.ndarray[double, ndim=2, mode='c'] u,
    np.ndarray[double, ndim=2, mode='c'] u_1,
    np.ndarray[double, ndim=2, mode='c'] u_2,
    np.ndarray[double, ndim=2, mode='c'] f,
    np.ndarray[double, ndim=1, mode='c'] x,
    np.ndarray[double, ndim=1, mode='c'] y,
    np.ndarray[double, ndim=1, mode='c'] t,
    double Cx2, double Cy2, double dt2):
    advance(&u[0,0], &u_1[0,0], &u_2[0,0], &f[0,0],
            &x[0], &y[0], &t[0],
            Cx2, Cy2, dt2,
            u.shape[0]-1, u.shape[1]-1, t.shape[0]-1)
    return u
