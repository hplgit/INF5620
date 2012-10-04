import numpy as np
cimport numpy as np
cimport cython
ctypedef np.float64_t DT    # data type

@cython.boundscheck(False)  # turn off array bounds check
@cython.wraparound(False)   # turn off negative indices (u[-1,-1])
cpdef advance(
    np.ndarray[DT, ndim=2, mode='c'] u,
    np.ndarray[DT, ndim=2, mode='c'] u_1,
    np.ndarray[DT, ndim=2, mode='c'] u_2,
    np.ndarray[DT, ndim=2, mode='c'] f,
    np.ndarray[DT, ndim=1, mode='c'] x,
    np.ndarray[DT, ndim=1, mode='c'] y,
    np.ndarray[DT, ndim=1, mode='c'] t,
    double Cx2, double Cy2, double dt2):

    cdef int Nx, Ny, i, j
    Nx = u.shape[0]-1
    Ny = u.shape[1]-1
    for i in xrange(1, Nx):
        for j in xrange(1, Ny):
            u[i,j] = 2*u_1[i,j] - u_2[i,j] + \
            Cx2*(u_1[i-1,j] - 2*u_1[i,j] + u_1[i+1,j]) + \
            Cy2*(u_1[i,j-1] - 2*u_1[i,j] + u_1[i,j+1]) + \
            dt2*f[i,j]
    # Boundary condition u=0
    j = 0
    for i in range(0, Nx+1): u[i,j] = 0
    j = Ny
    for i in range(0, Nx+1): u[i,j] = 0
    i = 0
    for j in range(0, Ny+1): u[i,j] = 0
    i = Nx
    for j in range(0, Ny+1): u[i,j] = 0
    return u
