#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         optSwTimes
# Date:         27.04.2024
# Author:       Dr. Pascal Schirmer
# Version:      V.1.0
# Copyright:    Pascal Schirmer
#######################################################################################################################
#######################################################################################################################

#######################################################################################################################
# Function Description
#######################################################################################################################
"""
This function calculates the optimal switching periods for OSM:
Inputs:     1) Mi:          modulation index (p.u.)
            2) fs:          switching frequency (Hz)
            3) K:           number of samples per period
            4) q:           pulse number
            5) N:           number of electrical periods
            6) dx:          switching times (normalised)
Outputs:    1) Tsk_c:       optimal continuous sample times
            2) Tsk:         optimal sample times
            3) tri:         optimal carrier for carrier based methods
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
# Function
#######################################################################################################################
def optSwTimes(Mi, fs, K, q, N, d0, d1, d2, d7):
    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Variables
    # ==============================================================================
    Tsk_c = np.zeros(len(d1))
    Tsk = np.zeros(q * N)
    ts = np.zeros(q * N)
    tri = np.zeros(len(d1))

    # ==============================================================================
    # Parameters
    # ==============================================================================
    c = np.pi / 3
    Ts = 1 / fs
    K_s = 0

    ###################################################################################################################
    # Pre-processing
    ###################################################################################################################
    # ==============================================================================
    # Stator Flux
    # ==============================================================================
    # ------------------------------------------
    # Coefficients
    # ------------------------------------------
    lmb_11 = (1 / 3 * (Mi / 4 * np.pi) ** 2 * (d0 ** 3 + d7 ** 3))
    lmb_12 = 1 / 3 * (c ** 2 + (Mi / 4 * np.pi) ** 2) * d2 ** 3 + (Mi / 4 * np.pi) ** 2 * d2 * d7 * (
                d7 + d2) - d2 ** 2 * c ** 2 * (d1 / 2 + d2) * (2 / 3 * d2 + d7)
    lmb_13 = 1 / 3 * c ** 2 * d1 ** 3 * (1 - d1) ** 2 + c ** 2 * d2 ** 2 * d0 * d1 * (
                d0 + d1) - c ** 2 * d0 * d1 ** 3 * (1 - d1 - d0) - 1 / 3 * c ** 2 * d1 ** 3 * d2 * (
                         1 - d1 - d2) + 1 / 2 * c ** 2 * d0 * d1 ** 2 * d2 * (2 * d0 - 1 + 2 * d1)

    # ------------------------------------------
    # Total Flux
    # ------------------------------------------
    lmb_v1 = (lmb_11 + lmb_12 + lmb_13) ** (1 / 3)
    lmb_sq = lmb_11 + lmb_12 + lmb_13

    # ------------------------------------------
    # Fix NaN
    # ------------------------------------------
    if np.sum(lmb_v1) == 0:
        lmb_v1 = np.ones(len(d1))
        lmb_sq = np.ones(len(d1))

    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    # ==============================================================================
    # Continuous
    # ==============================================================================
    for i in range(0, len(d1)):
        Tsk_c[i] = len(d1) * Ts / (lmb_v1[i] * np.sum(1 / lmb_v1))

    # ==============================================================================
    # Sampled
    # ==============================================================================
    temp = np.zeros(q*N)
    for i in range(0, q * N):
        temp[i] = np.sqrt(np.sum(lmb_sq[K*i:K*(i+1)])*1e-6)

    for i in range(0, q * N):
        Tsk[i] = (q*N*Ts) / (temp[i] * np.sum(1/temp)) / Ts
        ts[i] = K_s
        K_s = int(K * Tsk[i]) + K_s

    ###################################################################################################################
    # Post-processing
    ###################################################################################################################
    # ==============================================================================
    # Angle Correction
    # ==============================================================================
    Tsk_c = np.roll(Tsk_c, int(np.floor(0 / 360 * len(Tsk_c)))) / Ts

    # ==============================================================================
    # Carrier
    # ==============================================================================
    ende = 0
    for i in range(0, q * N):
        rise = np.linspace(1, -1, int(np.floor(K * Tsk[i]) / 2 + 1))
        fall = np.linspace(-1, 1, int(np.floor(K * Tsk[i]) / 2 + 2))
        pwm = np.concatenate((rise[0:-1], fall[0:-1]))
        if (ende + len(pwm)) < len(tri):
            tri[ende:ende + len(pwm)] = pwm
        else:
            tri[ende:len(tri)] = pwm[0:len(tri) - ende]
            break

        ende = ende + len(pwm)

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return [Tsk_c, Tsk, tri]
