#include <stdio.h>

#define idx(i,j) (i)*(Ny+1) + j
// (i) is essential, otherwise idx(i-1,j) becomes i-1*(Ny+1)+j!

void advance(double* u, double* u_1, double* u_2, double* f,
             double* x, double* y, double* t,
	     double Cx2, double Cy2, double dt2,
	     int Nx, int Ny, int N)
{
  int i, j;
  /* Scheme at interior points */
  for (i=1; i<=Nx-1; i++) {
    for (j=1; j<=Ny-1; j++) {
        u[idx(i,j)] = 2*u_1[idx(i,j)] - u_2[idx(i,j)] +
        Cx2*(u_1[idx(i-1,j)] - 2*u_1[idx(i,j)] + u_1[idx(i+1,j)]) +
        Cy2*(u_1[idx(i,j-1)] - 2*u_1[idx(i,j)] + u_1[idx(i,j+1)]) +
        dt2*f[idx(i,j)];
	if (i==1 && j==1) {
	  //printf("u[%d,%d]=%g, u_1[i,j]=%g %d\nu_1[i-1,j]=%g %d\nu_1[i+1,j]=%g %d\nu_1[i,j-1]=%g %d\nu_1[i,j+1]=%g %d\n", i,j,u[idx(i,j)], u_1[idx(i,j)], idx(i,j), u_1[idx(i-1,j)], idx(i-1,j), u_1[idx(i+1,j)], idx(i+1,j), u_1[idx(i,j-1)], idx(i,j-1), u_1[idx(i,j+1)], idx(i,j+1));
	}
    }
  }
  /* Boundary conditions */

  j = 0;  for (i=0; i<=Nx; i++) u[idx(i,j)] = 0;
  j = Ny; for (i=0; i<=Nx; i++) u[idx(i,j)] = 0;
  i = 0;  for (j=0; j<=Ny; j++) u[idx(i,j)] = 0;
  i = Nx; for (j=0; j<=Ny; j++) u[idx(i,j)] = 0;
}

