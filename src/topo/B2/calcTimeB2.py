#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         calcTimeB2
# Date:         14.08.2023
# Author:       Dr. Pascal A. Schirmer
# Version:      V.0.2
# Copyright:    Pascal Schirmer
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
import scipy.signal as sig


#######################################################################################################################
# Function
#######################################################################################################################
def calcTimeB2(t, s, e, Vdc, Mi, mdl, setupTopo, start, ende):
    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    timeAc = {}
    timeDc = {}

    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    # ==============================================================================
    # AC-Side
    # ==============================================================================
    # ------------------------------------------
    # Inverter Output
    # ------------------------------------------
    v_a0 = 0.5 * s * Vdc

    # ------------------------------------------
    # Filter Output
    # ------------------------------------------
    if setupTopo['outFilter'] == 0:
        v_L = v_a0
    else:
        _, v_L, _, = sig.lsim(mdl['SS']['Out'], v_a0, t, X0=v_a0[0])

    # ------------------------------------------
    # Load
    # ------------------------------------------
    # Voltage
    v_a = v_L - Mi * e

    # Current
    if setupTopo['wave'] == "con":
        _, i_a, _, = sig.lsim(mdl['SS']['Load'], v_a, t)
        i_a = i_a[start:ende]
    else:
        _, i_a, _, = sig.lsim(mdl['SS']['Load'], (v_a - np.mean(v_a)), t)
        i_a = i_a[start:ende]
        i_a = i_a - np.mean(i_a)

    # ==============================================================================
    # DC-Side
    # ==============================================================================
    # ------------------------------------------
    # Inverter Input
    # ------------------------------------------
    i_d_p = i_a * (1 + s[start:ende]) / 2
    i_d_m = i_a * (1 - s[start:ende]) / 2
    i_dc = i_d_p

    # ------------------------------------------
    # DC-Link
    # ------------------------------------------
    i_c = np.mean(i_d_p) - i_d_p
    _, v_dc, _, = sig.lsim(mdl['SS']['DC'], i_c, t[start:ende])

    # ------------------------------------------
    # Filter Input
    # ------------------------------------------
    if setupTopo['inpFilter'] == 0:
        v_in = v_dc
    else:
        _, v_in, _, = sig.lsim(mdl['SS']['Inp'], (v_dc - Vdc), t[start:ende])
        v_in = v_in + Vdc

    ###################################################################################################################
    # Post-Processing
    ###################################################################################################################
    # ==============================================================================
    # AC-Side
    # ==============================================================================
    timeAc['v_a0'] = v_a0[start:ende]
    timeAc['v_L'] = v_L[start:ende]
    timeAc['v_a'] = v_a[start:ende]
    timeAc['i_a'] = i_a

    # ==============================================================================
    # DC-Side
    # ==============================================================================
    timeDc['v_in'] = v_in
    timeDc['v_dc'] = v_dc
    timeDc['i_dc'] = i_dc
    timeDc['i_c'] = i_c
    timeDc['i_d_m'] = i_d_m
    timeDc['i_d_p'] = i_d_p

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return [timeAc, timeDc]
