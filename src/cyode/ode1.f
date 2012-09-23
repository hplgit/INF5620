      program ode2
      integer n_U0, n_t, n
      real*8 pi
      parameter (n_t=8000001)
      real*8 u(0:n_t-1), U0, dt, t(0:n_t-1)
      real*8 cpu0, cpu1
      external problem1, problem2, rk2
      U0 = 1.
      call cpu_time(cpu0)
      dt = 5./(n_t-1)
      t(0) = 0
      do n = 1, n_t-1
         t(n) = t(n-1) + dt
      end do
      call solver(problem1, U0, t, n_t, u, rk2)
      call cpu_time(cpu1)
      n = n_t-1
      write(*, 1000) 'u(', n, ')=', u(n)
      write(*, 2000) 'CPU time:', cpu1-cpu0
 1000 format(A, I8, A, F12.4)
 2000 format(A, F9.3)
      end


      subroutine solver(f, U0, t, n_t, u, method)
      integer n_u0, n_t
      real*8 U0, t(0:n_t-1), u(0:n_t-1)
      external f, method
      integer n
      real*8 dt 
      u(0) = U0

      do n = 0, n_t-2
         call method(u, n, t, f, n_t)
      end do
      return
      end

      subroutine rk2(u, n, t, f, n_t)
      integer n, n_t
      real*8 u(0:n_t-1), t(0:n_t-1)
      real*8 dt, un, dudt, K1, K2
      external f

      dt = t(n+1) - t(n)
      un = u(n)

      call f(dudt, un, t(n))

      K1 = dt*dudt
      un = u(n) + 0.5*K1

      call f(dudt, un, t(n) + 0.5*dt)

      K2 = dt*dudt
      u(n+1) = u(n) + K2
      return
      end

      subroutine problem1(dudt, u, t)
      dudt = -u + 1
      return
      end

      subroutine problem2(dudt, u, t)
      dudt = -u + exp(-2*t)
      return
      end
