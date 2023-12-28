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
from scipy.optimize import minimize, LinearConstraint, NonlinearConstraint


#######################################################################################################################
# Function
#######################################################################################################################
# Voltage amplitude function:
def ampl_sym_opp(alpha_a, k):
    '''
    Calkulates the voltage amplitudes of a quarter-wave symmetrical pulse-pattern for a three-phase system.
    Inputs:
    alpha_a: switching angles as numpy array
    k: highest Fourier-coefficient to be evaluated
    Output:
    u_ak(ki): Voltage amplitude over relative frequency ki
    i_ak(ki): Frequency in relation to fundamental frequency
    '''

    swi_alpha = [-2 if idx % 2 == 1 else 2 for idx,_ in enumerate(alpha_a)]
    i_ak = [x for x in range(1,k+1,2) if x % 3 != 0] # Fundamental and relevant harmonic components only
    u_ak = [0]*len(i_ak)

    for idxk, ikk in enumerate(i_ak):
        # Run over all relevant frequency components
        ka = np.cos(ikk*alpha_a)
        u_ak[idxk] = (np.dot(swi_alpha, ka)-1) / ikk

    return u_ak, i_ak

# Cost function:
def costfuntion(alpha):

    u_kc, i_kc = ampl_sym_opp(alpha,100)
    u_ak_sq = np.power(np.divide(u_kc[1:],i_kc[1:]), 2)
    wthd = np.sqrt(np.sum(u_ak_sq))
    cost = wthd

    return cost


def Opp(k_max, p0, Mi):
    '''
    Calculates the optimal switch- on and -off angles according to the following input parameters:
    k_max: Maximum order of harmonics, which are taken into account
    p0: Number of switching instances in one fundamental period
    Mi: Modulation index [0 ... 4/pi]
    '''

    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    p_deg = np.int16((p0 - 1)/2)        # Degrees of freedom (quarter wave symmetry)
    bmin = 0                            # Minimum angle between switching instances (tbd.)
    lbmin = bmin                        # Minium angle (lower bound for optimisation)
    ubmax = (np.pi - bmin)/2            # Maximum angle (upper bound for optimisation)
    alpha0 = [lbmin]*p_deg              # Starting point of alpha-values for optimisation
    ub = [ubmax]*p_deg                  # Upper bound list
    lb = [lbmin]*p_deg                  # Lower bound list
    bounds = list(zip(lb, ub))          # Bounds as a list of tuples

    ###################################################################################################################
    # Pre-processing
    ###################################################################################################################

    # Inequality constraints (linear):
    aineq = np.concatenate( (np.eye(p_deg-1, p_deg) -
                              np.concatenate((np.zeros((p_deg-1, 1)), np.eye(p_deg-1)), axis=1),
                              np.zeros((1, p_deg))), axis=0)
    bineq = np.append([-bmin]*(p_deg-1), [0])

    # Equality constraint:
    def eq_con(x_eq):

        u_fund, _ = ampl_sym_opp(x_eq, 1)
        return u_fund

    linear_constraint = LinearConstraint(aineq, lb=[-np.inf]*p_deg, ub=bineq)
    nonlinear_constraint = NonlinearConstraint(eq_con, lb=Mi, ub=Mi)

    ###################################################################################################################
    # Calculation
    ###################################################################################################################

    opt_result = minimize(costfuntion, alpha0, method='SLSQP', bounds=bounds,
                          constraints=[linear_constraint, nonlinear_constraint])


    ###################################################################################################################
    # Return
    ###################################################################################################################
    return opt_result

if __name__ == "__main__":

    u_k, i_k = ampl_sym_opp(np.array([np.pi/3]), 50)
    wthdc = costfuntion(np.array([0, 0, 0, 0]))

    result = Opp(100, 5, 0.5)
    print(result.message)
