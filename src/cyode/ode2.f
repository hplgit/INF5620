      program ode2
      integer n_U0, n_t, m, i, n, nperiods
C  9000001
      real*8 pi
      parameter (nperiods=30000)
      parameter (pi=3.14159, n_U0=2, m=n_U0-1, n_t=30*nperiods*2*pi-1)
      real*8 u(0:n_t-1, 0:m), U0(0:m), dt, t_end, t(0:n_t-1)
      real*8 work(4*n_U0), cpu0, cpu1
      external problem2, rk2
      t_end = 4*pi
      U0(0) = 0
      U0(1) = 1
      call cpu_time(cpu0)
      dt = t_end/(n_t-1)
      t(0) = 0
      do n = 1, n_t-1
         t(n) = t(n-1) + dt
      end do

      call solver(problem2, U0, n_U0, t, n_t, u, rk2, work)
      call cpu_time(cpu1)
C      do n = 0, n_t-1
C         do i = 0, m
C            write(*,*) 'u(', n, ',', i, ')=', u(n,i)
C         end do
C      end do
      n = n_t-1
      do i = 0, m
         write(*, 1000) 'u(', n, ',', i, ')=', u(n,i)
      end do
      write(*, 2000) 'CPU time:', cpu1-cpu0
 1000 format(A, I8, A, I2, A, F12.4)
 2000 format(A, F9.3)
      end


      subroutine solver(f, U0, n_U0, t, n_t, u, method, work)
      integer n_u0, n_t
      real*8 U0(0:n_U0-1), t(0:n_t-1), work(4*n_U0)
      real*8 u(0:n_t-1, 0:n_U0-1)
      external f, method

      integer n, m, i
      real*8 dt 
      m = n_U0-1

      do i = 0, m
         u(0,i) = U0(i)
      end do

      do n = 0, n_t-2
         call method(u, n, t, f, work(1), work(n_U0+1), work(2*n_U0+1), 
     &               work(3*n_U0+1), m, n_t)
      end do
      return
      end

      subroutine rk2(u, n, t, f, dudt, K1, K2, un, m, n_t)
      integer m, n, n_t
      real*8 u(0:n_t-1,0:m), t(0:n_t-1), dudt(0:m), K1(0:m), K2(0:m), 
     &       un(0:m)
      integer i
      real*8 dt
      external f

      dt = t(n+1) - t(n)
      do i = 0,m
         un(i) = u(n, i)
      end do

      call f(dudt, un, t(n), m)

      do i = 0,m
         K1(i) = dt*dudt(i)
      end do
      do i = 0,m
         un(i) = u(n, i) + 0.5*K1(i)
      end do

      call f(dudt, un, t(n) + 0.5*dt, m)

      do i = 0,m
         K2(i) = dt*dudt(i)
      end do
      do i = 0,m
         u(n+1,i) = u(n,i) + K2(i)
      end do
      return
      end

      subroutine problem2(dudt, u, t, m)
      integer m
      real*8 dudt(0:m), u(0:m), t
      dudt(0) = u(1)
      dudt(1) = -u(0)
      return
      end
