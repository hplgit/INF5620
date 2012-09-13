import scitools.std as plt
#import matplotlib.pyplot as plt
import sys
import odespy
import numpy as np

def f(u, t, w=1):
    # u is array of length 2 holding our [u, v]
    u, v = u
    return [v, -w**2*u]

def run_solvers_and_plot(solvers, timesteps_per_period=20,
                         num_periods=1, I=1, w=2*np.pi):
    P = 2*np.pi/w  # one period
    dt = P/timesteps_per_period
    N = num_periods*timesteps_per_period
    T = N*dt
    t_mesh = np.linspace(0, T, N+1)

    legends = []
    for solver in solvers:
        solver.set(f_kwargs={'w': w})
        solver.set_initial_condition([I, 0])
        u, t = solver.solve(t_mesh)

        # Make plots
        plt.figure(1)
        if len(t_mesh) <= 50:
            plt.plot(t, u[:,0])             # markers by default
        else:
            plt.plot(t, u[:,0], '-2')       # no markers
        plt.hold('on')
        legends.append(solver.__class__.__name__)
        plt.figure(2)
        if len(t_mesh) <= 50:
            plt.plot(u[:,0], u[:,1])        # markers by default
        else:
            plt.plot(u[:,0], u[:,1], '-2')  # no markers
        plt.hold('on')

    # Compare with exact solution plotted on a very fine mesh
    t_fine = np.linspace(0, T, 10001)
    u_e = I*np.cos(w*t_fine)
    v_e = -w*I*np.sin(w*t_fine)

    plt.figure(1)
    plt.plot(t_fine, u_e, '-') # avoid markers by spec. line type
    legends.append('exact')
    plt.legend(legends, loc='upper left')
    plt.xlabel('t');  plt.ylabel('u')
    plt.title('Time step: %g' % dt)
    plt.savefig('vb_%d_%d_u.png' % (timesteps_per_period, num_periods))
    plt.savefig('vb_%d_%d_u.pdf' % (timesteps_per_period, num_periods))
    plt.savefig('vb_%d_%d_u.eps' % (timesteps_per_period, num_periods))

    plt.figure(2)
    plt.plot(u_e, v_e, '-') # avoid markers by spec. line type
    plt.legend(legends, loc='lower right')
    plt.xlabel('u(t)');  plt.ylabel('v(t)')
    plt.title('Time step: %g' % dt)
    plt.savefig('vb_%d_%d_pp.png' % (timesteps_per_period, num_periods))
    plt.savefig('vb_%d_%d_pp.pdf' % (timesteps_per_period, num_periods))
    plt.savefig('vb_%d_%d_pp.eps' % (timesteps_per_period, num_periods))


solvers_theta = [
    odespy.ForwardEuler(f),
    # Implicit methods must use Newton solver to converge
    odespy.BackwardEuler   (f, nonlinear_solver='Newton'),
    odespy.MidpointImplicit(f, nonlinear_solver='Newton'),
    ]

solvers_RK = [odespy.RK2(f), odespy.RK4(f)]
solvers_CN = [odespy.MidpointImplicit(f, nonlinear_solver='Newton')]

if __name__ == '__main__':
    timesteps_per_period = 20
    solver_collection = 'theta'
    num_periods = 1
    try:
        timesteps_per_period = int(sys.argv[1])
        solver_collection = sys.argv[2]
        num_periods = int(sys.argv[3])
    except IndexError:
        pass # default values are ok
    run_solvers_and_plot(eval('solvers_' + solver_collection),
                         timesteps_per_period, num_periods)
    plt.show()
    raw_input()

