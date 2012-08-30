import nose.tools as nt
import decay_theta as decay
import numpy as np
import odespy

def test_odespy():
    """Compare odespy methods against decay.theta_rule."""
    I = 0.8; a = 1.2; T = 4; dt = 0.5
    prms0 = dict(eps_iter=1E-6, verbose=0)

    methods = [(0, odespy.ForwardEuler, {}), # (theta, classname, odespy prms)
               (0, odespy.ThetaRule, {'theta': 0})]
    for method in 'BackwardEuler', 'CrankNicolson', 0.5, 1:
        if method == 0.5:
            method = 'ThetaRule'
            theta = 0.5
        elif method == 1:
            method = 'ThetaRule'
            theta = 1
        elif method == 'CrankNicolson':
            theta = 0.5
        elif method == 'BackwardEuler':
            theta = 1

        for nlsolver in 'Picard', 'Newton':
            prms = prms0.copy()
            prms['theta'] = theta
            prms['nonlinear_solver'] = nlsolver
            methods.append((theta, eval('odespy.'+method), prms))
            if nlsolver == 'Newton':
                # add a version with Jacobian
                def myjac(u, t):
                    return -a
                #prms['jac'] = lambda u, t: -a
                prms2 = prms.copy()
                prms2['jac'] = myjac
                methods.append((theta, eval('odespy.'+method), prms2))
    for theta, odespy_solver, parameters in methods:
        u_d, t_d = decay.theta_rule(I, a, T, dt, theta)
        solver = odespy_solver(f=lambda u, t: -a*u, **parameters)
        solver.set_initial_condition(I)
        u, t = solver.solve(np.linspace(0, T, t_d.size))
        ok = np.allclose(u_d, u, rtol=10*prms0['eps_iter'], atol=prms0['eps_iter'])
        # str(solver) prints all parameters *different* from the defaults
        print str(solver), 'works' if ok else 'does not work'
        if not ok:
            print 'max deviation:', abs(u_d-u).max()
        assert ok
#test_odespy()
