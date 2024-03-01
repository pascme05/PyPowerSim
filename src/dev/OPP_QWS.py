#######################################################################################################################
# Import libs
#######################################################################################################################
# ==============================================================================
# Internal
# ==============================================================================

# ==============================================================================
# External
# ==============================================================================
import numpy as np
from scipy.optimize import LinearConstraint, NonlinearConstraint, basinhopping
import matplotlib.pyplot as plt


#######################################################################################################################
# Function
#######################################################################################################################
# Voltage amplitude function:
def ampl_sym_opp(alpha_a, kmax, sign):

    """
    Calculates the voltage amplitudes of a quarter-wave symmetrical pulse-pattern for a three-phase system.
    Inputs:
    alpha_a: switching angles as numpy array
    kmax: highest Fourier-coefficient to be evaluated
    Output:
    u_ak(ki): Voltage amplitude over relative frequency ki
    i_ak(ki): Frequency in relation to fundamental frequency
    """

    swi_alpha = [-2 if idx % 2 == 1 else 2 for idx, _ in enumerate(alpha_a)]
    i_ak = [x for x in range(1, kmax+1, 2) if x % 3 != 0]
    u_ak = [0]*len(i_ak)

    for idxk, ikk in enumerate(i_ak):
        # Run over all relevant frequency components
        ka = np.cos(ikk*alpha_a)
        u_ak[idxk] = sign * (np.dot(swi_alpha, ka)-1) / ikk

    return u_ak, i_ak


# Cost function:
def costfuntion(alpha, kmax):

    u_kc, i_kc = ampl_sym_opp(alpha, kmax, sign=1)
    u_ak_sq = np.power(np.divide(u_kc[1:], i_kc[1:]), 2)
    wthd = np.sqrt(np.sum(u_ak_sq))
    cost = wthd

    return cost


def print_fun(x, f, accepted):
    print("at minimum %.4f accepted %d" % (f, int(accepted)))


def Opp(kmax, p0, Mi, alpha0):
    """
    Calculates the optimal switch-on and -off angles according to the following input parameters:
    k_max: Maximum order of harmonics, which are taken into account
    p0: Number of switching instances in one fundamental period
    Mi: Modulation index [0 ... 4/pi]
    """

    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    p_deg = np.int16((p0 - 1)/2)        # Degrees of freedom (quarter wave symmetry)
    bmin = 0                            # Minimum angle between switching instances (tbd.)
    lbmin = bmin                        # Minimum angle (lower bound for optimisation)
    ubmax = (np.pi - bmin)/2            # Maximum angle (upper bound for optimisation)
    ub = [ubmax]*p_deg                  # Upper bound list
    lb = [lbmin]*p_deg                  # Lower bound list
    bounds = list(zip(lb, ub))          # Bounds as a list of tuples

    ###################################################################################################################
    # Pre-processing
    ###################################################################################################################
    # Inequality constraints (linear)
    aineq = np.concatenate((np.eye(p_deg-1, p_deg) -
                            np.concatenate((np.zeros((p_deg-1, 1)), np.eye(p_deg-1)), axis=1),
                            np.zeros((1, p_deg))), axis=0)
    bineq = np.append([-bmin]*(p_deg-1), [0])

    # Equality constraint
    def eq_con(x_eq):
        u_fund, _ = ampl_sym_opp(x_eq, 1, sign=1)
        return abs(u_fund[0])

    # Defining Solver Constraints
    linear_constraint = LinearConstraint(aineq, lb=[-np.inf]*p_deg, ub=bineq)
    nonlinear_constraint = NonlinearConstraint(eq_con, lb=Mi, ub=Mi)

    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    minimizer_kwargs = {"method": "SLSQP", "bounds": bounds, "constraints": [linear_constraint, nonlinear_constraint], "args": kmax}
    opt_result = basinhopping(costfuntion, alpha0, minimizer_kwargs=minimizer_kwargs, niter=200, callback=print_fun)

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return opt_result


if __name__ == "__main__":
    Num = 50
    p0 = 3
    kmax = 100
    thd = np.zeros((p0, Num))
    alpha = np.pi/3*np.ones((2*p0+1, p0, Num))
    Mi = np.linspace(0.0, 1, Num)
    for i in range(2, p0):
        p = 2*i + 1
        alpha0 = np.zeros(int((p - 1) / 2))
        for ii in range(0, Num):
            result = Opp(kmax, p, Mi[ii], alpha0)
            alpha[0:len(result.x), i, ii] = result.x
            alpha0 = result.x
            thd[i, ii] = result.fun
        plt.plot(thd[i, :])
        # plt.figure()
        # plt.plot(np.transpose(alpha[:, -1, :]))
    plt.show()
