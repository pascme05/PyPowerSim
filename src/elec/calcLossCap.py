#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         calcLossCap
# Date:         01.05.2024
# Author:       Dr. Pascal A. Schirmer
# Version:      V.1.0
# Copyright:    Pascal Schirmer
#######################################################################################################################
#######################################################################################################################

#######################################################################################################################
# Function Description
#######################################################################################################################
"""
This function calculates the losses of the capacitor devices.
Inputs:     1) t:       time vector (sec)
            2) i_c:     capacitor current (A)
            3) t_Tj:    core temperature of the capacitor (°C)
            4) para:    parameters of the switch
            5) setup:   all setup variables
Outputs:    1) out:     output array including capacitor
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
import pandas as pd
import numpy as np
from scipy import interpolate


#######################################################################################################################
# Function
#######################################################################################################################
def calcLossCap(t, i_c, Tj, para, setup):
    ###################################################################################################################
    # Initialisation
    ###################################################################################################################   
    # ==============================================================================
    # Variables
    # ==============================================================================
    dt = t[1] - t[0]
    f = np.linspace(0, int(1 / dt), int(len(i_c) / 2))
    fel = setup['Top']['fel']

    # ==============================================================================
    # Output
    # ==============================================================================
    out = pd.DataFrame(columns=['p_L'])

    ###################################################################################################################
    # Pre-Processing
    ###################################################################################################################
    # ==============================================================================
    # Get Fundamental Frequency
    # ==============================================================================
    idx = np.argmin(f - 2 * fel)

    # ==============================================================================
    # RMS Capacitor Current
    # ==============================================================================
    I_c_rms_sq = np.sum(i_c ** 2) / len(i_c)

    # ==============================================================================
    # Extract Parameters
    # ==============================================================================
    # ------------------------------------------
    # Constant
    # ------------------------------------------
    if setup['Par']['Elec']['CapMdl'] == "con" or setup['Par']['Elec']['CapMdl'] == "pwl":
        ESR = para['Cap']['Elec']['con']['ESR']

    # ------------------------------------------
    # Tabular
    # ------------------------------------------
    elif setup['Par']['Elec']['CapMdl'] == "tab":
        # Matrix 
        ESR_2d = interpolate.interp2d(para['Cap']['Elec']['vec']['Tj'].to_numpy(),
                                      para['Cap']['Elec']['vec']['f'].to_numpy(),
                                      para['Cap']['Elec']['tab']['ESR'].to_numpy(), kind='linear')

        # Fundamental Value
        ESR = ESR_2d(Tj, f[idx])

        # Static
        ESR = ESR[0]

    # ------------------------------------------
    # Default
    # ------------------------------------------
    else:
        ESR = para['Cap']['Elec']['con']['ESR']

    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    loss = ESR * I_c_rms_sq
    out['p_L'] = i_c ** 2 * (loss / (np.mean(i_c ** 2)))

    ###################################################################################################################
    # Post-Processing
    ###################################################################################################################
    # ==============================================================================
    # Scale Number Switches
    # ==============================================================================
    out['p_L'] = out['p_L'] * setup['Par']['Elec']['CapSeries'] * setup['Par']['Elec']['CapPara']

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return out
