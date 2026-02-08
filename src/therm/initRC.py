#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         initRC
# Date:         08.02.2026
# Author:       Dr. Pascal A. Schirmer
# Version:      V.1.0
# Copyright:    Pascal Schirmer
#######################################################################################################################
######################################################################################################################

#######################################################################################################################
# Function Description
#######################################################################################################################
"""
This function initialises the RC coefficient for thermal foster networks.
Inputs:     1) para:    all parameters used in the simulation
            2) setup:   includes all simulation variables
Outputs:    1) Rth:     thermal resistances (K/W)
            2) Cth:     thermal capacitance (Ws/K)
            The outputs are returned as a list containing:
            - Switch parameters (JA, DA, CA)
            - Capacitor parameters (JA)
            - Transformer parameters (PC, SC, CC, CA)
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
import pandas as pd


#######################################################################################################################
# Function
#######################################################################################################################
def initRC(para, setup):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("INFO: Initialise RC parameters")

    ###################################################################################################################
    # Loading Data
    ###################################################################################################################
    # ==============================================================================
    # Switches
    # ==============================================================================
    try:
        Rth_JC = para['Swi']['Ther']['vec']['Rth_JC'].values
        Cth_JC = para['Swi']['Ther']['vec']['Cth_JC'].values
        Rth_DC = para['Swi']['Ther']['vec']['Rth_DC'].values
        Cth_DC = para['Swi']['Ther']['vec']['Cth_DC'].values
        Rth_CA = para['Swi']['Ther']['vec']['Rth_CA'].values
        Cth_CA = para['Swi']['Ther']['vec']['Cth_CA'].values
    except KeyError:
        print("ERROR: Switch thermal parameters not found. Rth=0 K/W and Cth=1 Ws/K")
        Rth_JC, Cth_JC, Rth_DC, Cth_DC, Rth_CA, Cth_CA = np.array([0]), np.array([1]), np.array([0]), np.array([1]), np.array([0]), np.array([1])

    # ==============================================================================
    # Capacitor
    # ==============================================================================
    try:
        Rth_JC_cap = para['Cap']['Ther']['vec']['Rth_JC'].values
        Cth_JC_cap = para['Cap']['Ther']['vec']['Cth_JC'].values
        Rth_CA_cap = para['Cap']['Ther']['vec']['Rth_CA'].values
        Cth_CA_cap = para['Cap']['Ther']['vec']['Cth_CA'].values
    except (KeyError, TypeError):
        print("ERROR: Capacitor thermal parameters not found. Rth=0 K/W and Cth=1 Ws/K")
        Rth_JC_cap, Cth_JC_cap, Rth_CA_cap, Cth_CA_cap = np.array([0]), np.array([1]), np.array([0]), np.array([1])

    # ==============================================================================
    # Transformer
    # ==============================================================================
    try:
        Rth_PC_tra = para['Tra']['Ther']['vec']['Rth_PC'].values
        Cth_PC_tra = para['Tra']['Ther']['vec']['Cth_PC'].values
        Rth_SC_tra = para['Tra']['Ther']['vec']['Rth_SC'].values
        Cth_SC_tra = para['Tra']['Ther']['vec']['Cth_SC'].values
        Rth_CC_tra = para['Tra']['Ther']['vec']['Rth_CC'].values
        Cth_CC_tra = para['Tra']['Ther']['vec']['Cth_CC'].values
        Rth_CA_tra = para['Tra']['Ther']['vec']['Rth_CA'].values
        Cth_CA_tra = para['Tra']['Ther']['vec']['Cth_CA'].values

    except KeyError:
        print("ERROR: Transformer thermal parameters not found. Rth=0 K/W and Cth=1 Ws/K")
        Rth_PC_tra, Cth_PC_tra = np.array([0]), np.array([1])
        Rth_SC_tra, Cth_SC_tra = np.array([0]), np.array([1])
        Rth_CC_tra, Cth_CC_tra = np.array([0]), np.array([1])
        Rth_CA_tra, Cth_CA_tra = np.array([0]), np.array([1])

    ###################################################################################################################
    # Pre-Processing
    ###################################################################################################################
    # ==============================================================================
    # Switches
    # ==============================================================================
    Rth_JC = Rth_JC[np.logical_not(pd.isna(Rth_JC))]
    Cth_JC = Cth_JC[np.logical_not(pd.isna(Cth_JC))]
    Rth_DC = Rth_DC[np.logical_not(pd.isna(Rth_DC))]
    Cth_DC = Cth_DC[np.logical_not(pd.isna(Cth_DC))]
    Rth_CA = Rth_CA[np.logical_not(pd.isna(Rth_CA))]
    Cth_CA = Cth_CA[np.logical_not(pd.isna(Cth_CA))]

    # ==============================================================================
    # Capacitor
    # ==============================================================================
    Rth_JC_cap = Rth_JC_cap[np.logical_not(pd.isna(Rth_JC_cap))]
    Cth_JC_cap = Cth_JC_cap[np.logical_not(pd.isna(Cth_JC_cap))]
    Rth_CA_cap = Rth_CA_cap[np.logical_not(pd.isna(Rth_CA_cap))]
    Cth_CA_cap = Cth_CA_cap[np.logical_not(pd.isna(Cth_CA_cap))]

    # ==============================================================================
    # Transformer
    # ==============================================================================
    Rth_PC_tra = Rth_PC_tra[np.logical_not(pd.isna(Rth_PC_tra))]
    Cth_PC_tra = Cth_PC_tra[np.logical_not(pd.isna(Cth_PC_tra))]
    Rth_SC_tra = Rth_SC_tra[np.logical_not(pd.isna(Rth_SC_tra))]
    Cth_SC_tra = Cth_SC_tra[np.logical_not(pd.isna(Cth_SC_tra))]
    Rth_CC_tra = Rth_CC_tra[np.logical_not(pd.isna(Rth_CC_tra))]
    Cth_CC_tra = Cth_CC_tra[np.logical_not(pd.isna(Cth_CC_tra))]
    Rth_CA_tra = Rth_CA_tra[np.logical_not(pd.isna(Rth_CA_tra))]
    Cth_CA_tra = Cth_CA_tra[np.logical_not(pd.isna(Cth_CA_tra))]

    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    if setup['Par']['Ther']['Heatsink'] == 1 and setup['Par']['Ther']['Coupling'] == 0:
        # ==============================================================================
        # Switches
        # =============================================================================
        Rth_JA = np.concatenate((Rth_JC, Rth_CA))
        Cth_JA = np.concatenate((Cth_JC, Cth_CA))
        Rth_DA = np.concatenate((Rth_DC, Rth_CA))
        Cth_DA = np.concatenate((Cth_DC, Cth_CA))

        # ==============================================================================
        # Capacitor
        # =============================================================================
        Rth_JA_cap = np.concatenate((Rth_JC_cap, Rth_CA_cap))
        Cth_JA_cap = np.concatenate((Cth_JC_cap, Cth_CA_cap))

        # ==============================================================================
        # Transformer
        # =============================================================================
        Rth_PA_tra = np.concatenate((Rth_PC_tra, Rth_CA_tra))
        Cth_PA_tra = np.concatenate((Cth_PC_tra, Cth_CA_tra))
        Rth_SA_tra = np.concatenate((Rth_SC_tra, Rth_CA_tra))
        Cth_SA_tra = np.concatenate((Cth_SC_tra, Cth_CA_tra))
        Rth_CA_tra = np.concatenate((Rth_CC_tra, Rth_CA_tra))
        Cth_CA_tra = np.concatenate((Cth_CC_tra, Cth_CA_tra))
    else:
        Rth_JA = Rth_JC
        Cth_JA = Cth_JC
        Rth_DA = Rth_DC
        Cth_DA = Cth_DC
        Rth_JA_cap = Rth_JC_cap
        Cth_JA_cap = Cth_JC_cap
        Rth_PA_tra = Rth_PC_tra
        Cth_PA_tra = Cth_PC_tra
        Rth_SA_tra = Rth_SC_tra
        Cth_SA_tra = Cth_SC_tra
        Rth_CA_tra = Rth_CC_tra
        Cth_CA_tra = Cth_CC_tra
        
    ###################################################################################################################
    # Return
    ###################################################################################################################
    return [Rth_JA, Cth_JA, Rth_DA, Cth_DA, Rth_CA, Cth_CA, Rth_JA_cap, Cth_JA_cap,
            Rth_PA_tra, Cth_PA_tra, Rth_SA_tra, Cth_SA_tra, Rth_CA_tra, Cth_CA_tra]
