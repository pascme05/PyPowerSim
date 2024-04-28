#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         oppPWM
# Date:         27.04.2024
# Author:       Dr. Daniel Glose
# Version:      V.1.0
# Copyright:    Pascal Schirmer, Daniel Glose
#######################################################################################################################
#######################################################################################################################

#######################################################################################################################
# Function Description
#######################################################################################################################
"""
This function calculates the optimal switch-on and -off angles according to the following input parameters:
Inputs:     1) kmax:        maximum order of harmonics, which are taken into account
            2) p0:          pulse-number
            3) Mi:          modulation index (p.u.)
            4) sym:         quarter or half-wave symmetry
            5) setup:       includes all simulation variables
Outputs:    1) ang_total:   total vector of switching angles (rad)
            2) val_total:   total vector of switching states (-)
            3) wthd:        weighted total harmonic distortion
"""

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


#######################################################################################################################
# Additional Functions
#######################################################################################################################
# ==============================================================================
# Voltage amplitude function
# ==============================================================================
# ------------------------------------------
# B6
# ------------------------------------------
def ampl_sym_B6(alpha_a, kmax):
    # Switching angles
    swi_alpha = [-2 if idx % 2 == 1 else 2 for idx, _ in enumerate(alpha_a)]

    # Harmonics
    i_ak = [x for x in range(1, kmax + 1, 2) if x % 3 != 0]

    # Init
    u_ak = [0]*len(i_ak)

    # Relevant frequency components
    for idxk, ikk in enumerate(i_ak):
        ka = np.cos(ikk*alpha_a)
        u_ak[idxk] = (np.dot(swi_alpha, ka)-1) / ikk

    return u_ak, i_ak


# ------------------------------------------
# B4
# ------------------------------------------
def ampl_sym_B4(alpha_a, kmax):
    # Switching angles
    swi_alpha = [-1 if idx % 2 == 1 else 1 for idx, _ in enumerate(alpha_a)]

    # Harmonics
    i_ak = [x for x in range(1, kmax + 1, 2)]

    # Init
    u_ak = [0]*len(i_ak)

    # Relevant frequency components
    for idxk, ikk in enumerate(i_ak):
        ka = np.cos(ikk*alpha_a)
        u_ak[idxk] = (np.dot(swi_alpha, ka)) / ikk

    return u_ak, i_ak


# ------------------------------------------
# B2
# ------------------------------------------
def ampl_sym_B2(alpha_a, kmax):
    # Switching angles
    swi_alpha = [-2 if idx % 2 == 1 else 2 for idx, _ in enumerate(alpha_a)]

    # Harmonics
    i_ak = [x for x in range(1, kmax + 1, 2)]

    # Init
    u_ak = [0]*len(i_ak)

    # Relevant frequency components
    for idxk, ikk in enumerate(i_ak):
        ka = np.cos(ikk*alpha_a)
        u_ak[idxk] = (np.dot(swi_alpha, ka)-1) / ikk

    return u_ak, i_ak


# ==============================================================================
# Cost function:
# ==============================================================================
def costfunction(alpha, kmax, topo):
    # ------------------------------------------
    # B2
    # ------------------------------------------
    if topo == 'B2':
        u_kc, i_kc = ampl_sym_B2(alpha, kmax)

    # ------------------------------------------
    # B4
    # ------------------------------------------
    elif topo == 'B4':
        u_kc, i_kc = ampl_sym_B4(alpha, kmax)

    # ------------------------------------------
    # B6
    # ------------------------------------------
    else:
        u_kc, i_kc = ampl_sym_B6(alpha, kmax)

    # ------------------------------------------
    # WTHD
    # ------------------------------------------
    u_ak_sq = np.power(np.divide(u_kc[1:], i_kc[1:]), 2)
    wthd = np.sqrt(np.sum(u_ak_sq))

    return wthd


# ==============================================================================
# Equality constraint
# ==============================================================================
# ------------------------------------------
# B6
# ------------------------------------------
def eq_B6(x_eq):
    swi_alpha = [-2 if idx % 2 == 1 else 2 for idx, _ in enumerate(x_eq)]
    ka = np.cos(x_eq)
    u1 = abs(np.sum((np.dot(swi_alpha, ka) - 1)))

    return u1


# ------------------------------------------
# B6
# ------------------------------------------
def eq_B4(x_eq):
    swi_alpha = [-1 if idx % 2 == 1 else 1 for idx, _ in enumerate(x_eq)]
    u1 = np.sum(np.dot(swi_alpha, np.cos(x_eq)))

    return u1


# ------------------------------------------
# B2
# ------------------------------------------
def eq_B2(x_eq):
    swi_alpha = [-2 if idx % 2 == 1 else 2 for idx, _ in enumerate(x_eq)]
    ka = np.cos(x_eq)
    u1 = abs(np.sum((np.dot(swi_alpha, ka) - 1)))

    return u1


#######################################################################################################################
# Function
#######################################################################################################################
def oppPWM(kmax, p0, Mi, sym, setupTopo):
    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    p_deg = np.int16((p0 - 1)/2)                                                                                        # Degrees of freedom (quarter wave symmetry)
    bmin = 0                                                                                                            # Minimum angle between switching instances (tbd.)
    lbmin = bmin                                                                                                        # Minimum angle (lower bound for optimisation)
    ubmax = (np.pi - bmin)/2                                                                                            # Maximum angle (upper bound for optimisation)
    alpha0 = [lbmin]*p_deg                                                                                              # Starting point of alpha-values for optimisation
    ub = [ubmax]*p_deg                                                                                                  # Upper bound list
    lb = [lbmin]*p_deg                                                                                                  # Lower bound list
    Mi_iter = np.linspace(0, Mi, 20)                                                                                    # Iteration for Modulation Index
    bounds = list(zip(lb, ub))                                                                                          # Bounds as a list of tuples
    opt_result = []

    ###################################################################################################################
    # Pre-processing
    ###################################################################################################################
    # ==============================================================================
    # Inequality constraints
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
        # ==============================================================================
        # Constraint
        # ==============================================================================
        # ------------------------------------------
        # B2
        # ------------------------------------------
        if setupTopo['sourceType'] == 'B2':
            nonlinear_constraint = NonlinearConstraint(eq_B2, lb=Mi, ub=Mi)

        # ------------------------------------------
        # B2
        # ------------------------------------------
        elif setupTopo['sourceType'] == 'B4':
            nonlinear_constraint = NonlinearConstraint(eq_B4, lb=Mi, ub=Mi)

        # ------------------------------------------
        # B2
        # ------------------------------------------
        else:
            nonlinear_constraint = NonlinearConstraint(eq_B6, lb=Mi, ub=Mi)

        # ==============================================================================
        # Solve
        # ==============================================================================
        minimizer_kwargs = {"method": "SLSQP", "bounds": bounds, "args": (kmax, setupTopo['sourceType']),
                            "constraints": [linear_constraint, nonlinear_constraint]}
        opt_result = basinhopping(costfunction, alpha0, minimizer_kwargs=minimizer_kwargs, niter=100)

        # ==============================================================================
        # Out
        # ==============================================================================
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
