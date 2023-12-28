#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         oppPWM
# Date:         13.11.2023
# Author:       Dr. Daniel Glose
# Version:      V.0.2
# Copyright:    Pascal Schirmer, Daniel Glose
#######################################################################################################################
#######################################################################################################################

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
from scipy.optimize import minimize, LinearConstraint, NonlinearConstraint, basinhopping


#######################################################################################################################
# Additional Functions
#######################################################################################################################
# ==============================================================================
# Voltage amplitude function
# ==============================================================================
def ampl_sym_opp(alpha_a, k):
    # ------------------------------------------
    # Theory
    # ------------------------------------------
    """
    Calculates the voltage amplitudes of a quarter-wave symmetrical pulse-pattern for a three-phase system.
    Inputs:
    alpha_a: switching angles as numpy array
    k: highest Fourier-coefficient to be evaluated
    Output:
    u_ak(ki): Voltage amplitude over relative frequency ki
    i_ak(ki): Frequency in relation to fundamental frequency
    """

    # ------------------------------------------
    # Theory
    # ------------------------------------------
    swi_alpha = [-2 if idx % 2 == 1 else 2 for idx, _ in enumerate(alpha_a)]
    #i_ak = [x for x in range(1, k+1, 2)]
    i_ak = [x for x in range(1,k+1,2) if x % 3 != 0]
    u_ak = [0]*len(i_ak)

    # ------------------------------------------
    # Relevant frequency components
    # ------------------------------------------
    for idxk, ikk in enumerate(i_ak):
        ka = np.cos(ikk*alpha_a)
        u_ak[idxk] = (np.dot(swi_alpha, ka)-1) / ikk

    return u_ak, i_ak


# ==============================================================================
# Cost function:
# ==============================================================================
def costfuntion(alpha):

    u_kc, i_kc = ampl_sym_opp(alpha, 100)
    u_ak_sq = np.power(np.divide(u_kc[1:], i_kc[1:]), 2)
    wthd = np.sqrt(np.sum(u_ak_sq))

    return wthd


# ==============================================================================
# Equality constraint
# ==============================================================================
def eq_con(x_eq):

    u_fund, _ = ampl_sym_opp(x_eq, 1)

    return u_fund[0]


#######################################################################################################################
# Function
#######################################################################################################################
def oppPWM(k_max, p0, Mi, sym):
    ###################################################################################################################
    # Theory
    ###################################################################################################################
    """
    Calculates the optimal switch-on and -off angles according to the following input parameters:
    k_max: Maximum order of harmonics, which are taken into account
    p0: Number of switching instances in one fundamental period
    Mi: Modulation index [0 ... 4/pi]
    sym: symmetrical boundary, e.g. quarter-wave
    """

    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    p_deg = np.int16((p0 - 1)/2)        # Degrees of freedom (quarter wave symmetry)
    bmin = 0                            # Minimum angle between switching instances (tbd.)
    lbmin = bmin                        # Minimum angle (lower bound for optimisation)
    ubmax = (np.pi - bmin)/2            # Maximum angle (upper bound for optimisation)
    alpha0 = [lbmin]*p_deg              # Starting point of alpha-values for optimisation
    ub = [ubmax]*p_deg                  # Upper bound list
    lb = [lbmin]*p_deg                  # Lower bound list
    Mi_iter = np.linspace(0, Mi, 20)
    bounds = list(zip(lb, ub))          # Bounds as a list of tuples

    ###################################################################################################################
    # Pre-processing
    ###################################################################################################################
    # ==============================================================================
    # Init
    # ==============================================================================
    aineq = np.concatenate((np.eye(p_deg-1, p_deg) -
                            np.concatenate((np.zeros((p_deg-1, 1)), np.eye(p_deg-1)), axis=1),
                            np.zeros((1, p_deg))), axis=0)
    bineq = np.append([-bmin]*(p_deg-1), [0])

    # ==============================================================================
    # Equality constraints
    # ==============================================================================
    linear_constraint = LinearConstraint(aineq, lb=[-np.inf]*p_deg, ub=bineq)

    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    for i in range(0, len(Mi_iter)):
        nonlinear_constraint = NonlinearConstraint(eq_con, lb=Mi_iter[i], ub=Mi_iter[i])
        minimizer_kwargs = {"method": "SLSQP", "bounds": bounds,
                            "constraints": [linear_constraint, nonlinear_constraint]}

        opt_result = basinhopping(costfuntion, alpha0, minimizer_kwargs=minimizer_kwargs, niter=20)
        #opt_result = minimize(costfuntion, alpha0, method='SLSQP', bounds=bounds,
         #                     constraints=[linear_constraint, nonlinear_constraint])
        alpha0 = opt_result.x

    ###################################################################################################################
    # Post-processing
    ###################################################################################################################
    # ==============================================================================
    # Angles
    # ==============================================================================
    # ------------------------------------------
    # Quarter Wave
    # ------------------------------------------
    if sym == 4:
        ang = opt_result.x
        ang_total = np.concatenate((ang, np.pi - np.flip(ang), np.pi + ang, 2 * np.pi - np.flip(ang)), axis=0)

    # ------------------------------------------
    # Default
    # ------------------------------------------
    else:
        ang = opt_result.x
        ang_total = np.concatenate((ang, np.pi - np.flip(ang), np.pi + ang, 2 * np.pi - np.flip(ang)), axis=0)

    # ==============================================================================
    # Direction
    # ==============================================================================
    # ------------------------------------------
    # Quarter Wave
    # ------------------------------------------
    if sym == 4:
        val = np.ones(np.size(ang))
        val_total = np.concatenate((val, val, -val, -val), axis=0)

    # ------------------------------------------
    # Default
    # ------------------------------------------
    else:
        val = np.ones(np.size(ang))
        val_total = np.concatenate((val, val, -val, -val), axis=0)

    # ==============================================================================
    # Distortion
    # ==============================================================================
    wthd = opt_result.fun

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return [ang_total, val_total, wthd]
