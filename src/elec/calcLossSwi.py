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
This function calculates the losses of the switching devices.
Inputs:     1) i_G:     gate signal
            2) i_T:     transistor current (A)
            3) i_D:     diode current (A)
            4) v_T:     transistor voltage (V)
            5) v_D:     diode voltage (V)
            6) t_Tj:    junction temperature of the switch (°C)
            7) para:    parameters of the switch
            8) setup:   all setup variables
Outputs:    1) out:     output array including transistor and diode quantities
"""

#######################################################################################################################
# Import libs
#######################################################################################################################
# ==============================================================================
# Internal
# ==============================================================================
from src.general.helpFnc import zoh_easy

# ==============================================================================
# External
# ==============================================================================
import numpy as np
import pandas as pd


#######################################################################################################################
# Function
#######################################################################################################################
def calcLossSwi(i_G, i_T, i_D, v_T, v_D, t_Tj, para, setup):
    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Parameters
    # ==============================================================================
    fs = setup['Par']['PWM']['fs']

    # ==============================================================================
    # Output
    # ==============================================================================
    out = pd.DataFrame(columns=['p_T_c', 'p_T_s', 'p_T', 'p_D_c', 'p_D_s', 'p_D', 'p_L'])

    ###################################################################################################################
    # Pre-Processing
    ###################################################################################################################
    # ==============================================================================
    # Extract Parameters
    # ==============================================================================
    # ------------------------------------------
    # Constant
    # ------------------------------------------
    if setup['Par']['Elec']['SwiMdl'] == "con":
        # IGBT
        if setup['Par']['Elec']['SwiType'] == "IGBT":
            Eon = para['Swi']['Elec']['con']['Eon'] * np.ones(np.size(i_T))
            Eoff = para['Swi']['Elec']['con']['Eoff'] * np.ones(np.size(i_T))
            Erec = para['Swi']['Elec']['con']['Erec'] * np.ones(np.size(i_T))

        # MOSFET
        else:
            Eon = 0.50 * np.max(np.abs(v_T)) * np.max(np.abs(i_T)) * (para['Swi']['Elec']['con']['tr'] + para['Swi']['Elec']['con']['tf']) + 0.5 * para['Swi']['Elec']['con']['Coss'] * np.max(np.abs(v_T)) ** 2 * np.ones(np.size(i_T))
            Eoff = 0.50 * np.max(np.abs(v_T)) * np.max(np.abs(i_T)) * (para['Swi']['Elec']['con']['tr'] + para['Swi']['Elec']['con']['tf']) + 0.5 * para['Swi']['Elec']['con']['Coss'] * np.max(np.abs(v_T)) ** 2 * np.ones(np.size(i_T))
            Erec = 0.25 * para['Swi']['Elec']['con']['Qrr'] * np.max(np.abs(v_D)) * np.ones(np.size(i_T)) + 0.5 * para['Swi']['Elec']['con']['Crss'] * np.max(np.abs(v_D)) ** 2 * np.ones(np.size(i_T))

    # ------------------------------------------
    # Piece-wise-linear
    # ------------------------------------------
    elif setup['Par']['Elec']['SwiMdl'] == "pwl":
        # IGBT
        if setup['Par']['Elec']['SwiType'] == "IGBT":
            Eon = para['Swi']['Elec']['con']['Eon'] * np.abs(i_T) / para['Swi']['Elec']['con']['Inom'] * np.max(np.abs(v_T)) / para['Swi']['Elec']['con']['Vnom']
            Eoff = para['Swi']['Elec']['con']['Eoff'] * np.abs(i_T) / para['Swi']['Elec']['con']['Inom'] * np.max(np.abs(v_T)) / para['Swi']['Elec']['con']['Vnom']
            Erec = para['Swi']['Elec']['con']['Erec'] * np.abs(i_D) / para['Swi']['Elec']['con']['Inom'] * np.max(np.abs(v_D)) / para['Swi']['Elec']['con']['Vnom']

        # MOSFET
        else:
            Eon = 0.50 * np.max(np.abs(v_T)) * np.abs(i_T) * (para['Swi']['Elec']['con']['tr'] + para['Swi']['Elec']['con']['tf']) + para['Swi']['Elec']['con']['Qrr'] * np.max(np.abs(v_T))
            Eoff = 0.50 * np.max(np.abs(v_T)) * np.abs(i_T) * (para['Swi']['Elec']['con']['tr'] + para['Swi']['Elec']['con']['tf'])
            Erec = 0.25 * para['Swi']['Elec']['con']['Qrr'] * np.max(np.abs(v_D)) * np.ones(np.size(i_T))

    # ------------------------------------------
    # Tabular
    # ------------------------------------------
    else:
        # IGBT
        if setup['Par']['Elec']['SwiType'] == "IGBT" or setup['Par']['PWM']['swloss'] == 0:
            Eon = para['Swi']['Elec']['tab']['Eon_2d']((t_Tj, abs(i_T))) * np.max(np.abs(v_T)) / para['Swi']['Elec']['con']['Vnom']
            Eoff = para['Swi']['Elec']['tab']['Eoff_2d']((t_Tj, abs(i_T))) * np.max(np.abs(v_T)) / para['Swi']['Elec']['con']['Vnom']
            Erec = para['Swi']['Elec']['tab']['Erec_2d']((t_Tj, abs(i_D))) * np.max(np.abs(v_D)) / para['Swi']['Elec']['con']['Vnom']

        # MOSFET
        else:
            IGon = (para['Swi']['Elec']['con']['Vg'] - para['Swi']['Elec']['con']['Vpl']) / para['Swi']['Elec']['con']['Rg']
            IGoff = para['Swi']['Elec']['con']['Vpl'] / para['Swi']['Elec']['con']['Rg']

            Eon = para['Swi']['Elec']['tab']['Eon_2d']((t_Tj, abs(v_T)))
            Eoff = para['Swi']['Elec']['tab']['Eon_2d']((t_Tj, abs(v_T)))
            Erec = para['Swi']['Elec']['tab']['Eon_2d']((t_Tj, abs(v_D)))
            Crss1 = para['Swi']['Elec']['tab']['Crss_2d']((t_Tj, abs(v_T)))
            Crss2 = para['Swi']['Elec']['tab']['Crss_2d']((t_Tj, np.max(np.abs(v_T))*np.ones(len(v_T))))

            tf1 = (np.max(np.abs(v_T)) - np.abs(v_T)) * (Crss1 / IGon)
            tf2 = (np.max(np.abs(v_T)) - np.abs(v_T)) * (Crss2 / IGon)
            tr1 = (np.max(np.abs(v_T)) - np.abs(v_T)) * (Crss1 / IGoff)
            tr2 = (np.max(np.abs(v_T)) - np.abs(v_T)) * (Crss2 / IGoff)
            tf = (tf1 + tf2) / 2
            tr = (tr1 + tr2) / 2

            Eon = Eon + 0.5 * np.abs(i_T) * np.max(abs(v_T)) * (tr + tf)
            Eoff = Eoff + 0.5 * np.abs(i_T) * np.max(abs(v_T)) * (tr + tf)
            Erec = Erec + 0.25 * para['Swi']['Elec']['con']['Qrr'] * np.max(np.abs(v_D))

    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    # ==============================================================================
    # Transistor
    # ==============================================================================
    # ------------------------------------------
    # Conduction
    # ------------------------------------------
    out['p_T_c'] = i_T * v_T

    # ------------------------------------------
    # Switching
    # ------------------------------------------
    out['p_T_s'] = (zoh_easy(Eon, np.diff(i_G, prepend=0)) + np.roll(zoh_easy(Eoff, np.diff(i_G, append=0) * (-1)), 1)) * fs * 2

    # ==============================================================================
    # Diode
    # ==============================================================================
    # ------------------------------------------
    # Conduction
    # ------------------------------------------
    out['p_D_c'] = i_D * v_D

    # ------------------------------------------
    # Switching
    # ------------------------------------------
    out['p_D_s'] = (zoh_easy(Erec, np.diff(i_G, prepend=0)) + np.roll(zoh_easy(Erec, np.diff(i_G, append=0) * (-1)), 1)) * fs

    # ==============================================================================
    # Total
    # ==============================================================================
    out['p_T'] = out['p_T_c'] + out['p_T_s']
    out['p_D'] = out['p_D_c'] + out['p_D_s']
    out['p_L'] = out['p_T'] + out['p_D']

    ###################################################################################################################
    # Post-Processing
    ###################################################################################################################
    # ==============================================================================
    # Scale Number Switches
    # ==============================================================================
    out['p_L'] = out['p_L'] * setup['Par']['Elec']['SwiSeries'] * setup['Par']['Elec']['SwiPara']

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return out
