      subroutine solver(I_func, V_func, f_func, c, Nx, Ny, dt,
     &                  dt_safety_factor,
     &                  x, y, t, N, u, u_1, u_2, f)
      integer Nx, Ny, N
      real*8 c, dt, dt_safety_factor
      real*8 x(0:Nx), y(0:Ny), t(0:N)
      real*8 u(0:Nx, 0:Ny), u_1(0:Nx, 0:Ny), u_2(0:Nx, 0:Ny),
     &       f(0:Nx, 0:Ny)
      real*8 I_func, V_func, f_func
      external I_func, V_func, f_func
      real*8 dx, dy, stability_limit, Cx2, Cy2, dt2
      integer i, j, n

C     Make coordinate vectors
      dx = x(1) - x(0)
      dy = y(1) - y(0)
      stability_limit = 1/c*(1/sqrt(1/dx**2 + 1/dy**2))
      if dt .le. 0 then
          dt = dt_safety_factor*stability_limit
      Cx2 = (c*dt/dx)**2
      Cy2 = (c*dt/dy)**2
      dt2 = dt**2

C     Compute initial condition
      do j = 0, Ny
         do i = 0, Nx
            u_1(i, j) = I_func(x(i), y(j))
         end do
      end do

C     Special formula for the first time step
      n = 0
      do j = 0, Ny
         do i = 0, Nx
            u(i, j) = u_1(i, j) + dt*V_func(x(i), y(j)) +
     &                0.5*Cx2*(u_1(i-1,j) - 2*u_1(i,j) + u_1(i+1,j)) +
     &                0.5*Cy2*(u_1(i,j-1) - 2*u_1(i,j) + u_1(i,j+1)) +
     &                0.5*dt2*f_func(x[i], y[j], t[n])
     &
         end do
      end do

C     Ensure u=0 on the boundary
      call bc_u0(u, Nx, Ny)

C     Switch arrays before next time step
      do j = 0, Ny
         do i = 0, Nx
            u_2(i, j) = u_1(i, j)
            u_1(i, j) = u(i, j)
         end do
      end do

C     Time loop
      do n = 1, N-1
C        Fill the f array with values
         do j = 0, Ny
            do i = 0, Nx
               f(i, j) = f_func(x(i), y(j), t(n))
            end do
         end do

         call advance(u, u_1, u_2, f, Cx2, Cy2, dt2, Nx, Ny)

         call write_to_file(u, Nx, Ny, n, t(n))

C        Switch arrays before next time step
         do j = 0, Ny
            do i = 0, Nx
               u_2(i, j) = u_1(i, j)
               u_1(i, j) = u(i, j)
            end do
         end do

      subroutine bc_u0(u, Nx, Ny)
C     Set u=0 on the boundary
      integer Nx, Ny, i, j
      real*8 u(0:Nx, 0:Ny)
      j = 0
      do i = 0, Nx
         u(i, j) = 0
      end do
      j = Ny
      do i = 0, Nx
         u(i, j) = 0
      end do
      i = 0
      do j = 0, Ny
         u(i, j) = 0
      end do
      i = Nx
      do j = 0, Ny
         u(i, j) = 0
      end do
      return
      end

      subroutine write_to_file(u, Nx, Ny, n, time)
      integer Nx, Ny, n
      real*8 time, u(0:Nx, 0:Ny)



