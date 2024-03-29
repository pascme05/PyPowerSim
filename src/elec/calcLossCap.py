#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         calcLossCap
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
import pandas as pd
import numpy as np
from scipy import interpolate


#######################################################################################################################
# Function
#######################################################################################################################
def calcLossCap(t, i_c, Tj, para, setupPara, setupTopo):
    ###################################################################################################################
    # Initialisation
    ###################################################################################################################   
    # ==============================================================================
    # Variables
    # ==============================================================================
    dt = t[1] - t[0]
    f = np.linspace(0, int(1 / dt), int(len(i_c) / 2))
    fel = setupTopo['fel']

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
    if setupPara['Elec']['CapMdl'] == "con" or setupPara['Elec']['CapMdl'] == "pwl":
        ESR = para['Cap']['Elec']['con']['ESR']

    # ------------------------------------------
    # Tabular
    # ------------------------------------------
    elif setupPara['Elec']['CapMdl'] == "tab":
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
    out['p_L'] = out['p_L'] * setupPara['Elec']['CapSeries'] * setupPara['Elec']['CapPara']

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return out
