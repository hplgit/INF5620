      subroutine advance(u, u_1, u_2, f, x, y, t, Cx2, Cy2, dt2,
     &                   Nx, Ny, N)
Cf2py intent(c) advance
      integer Nx, Ny, N
      real*8 u(0:Nx,0:Ny), u_1(0:Nx,0:Ny), u_2(0:Nx,0:Ny)
      real*8 f(0:Nx, 0:Ny), x(0:Nx), y(0:Ny), t(0:N), Cx2, Cy2, dt2
Cf2py intent(in, out) u
Cf2py intent(c) u, u_1, u_2, f, x, y, t, Cx2, Cy2, dt2, Nx, Ny, N
      return
      end


