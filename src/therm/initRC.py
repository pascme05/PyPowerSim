#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         initRC
# Date:         01.04.2023
# Author:       Dr. Pascal A. Schirmer
# Version:      V.0.1
# Copyright:    Pascal Schirmer
#######################################################################################################################
######################################################################################################################

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
def initRC(para, setupPara):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("INFO: Initialise RC parameters")

    ###################################################################################################################
    # Loading Data
    ###################################################################################################################
    Rth_JC = para['Swi']['Ther']['vec']['Rth_JC'].values
    Cth_JC = para['Swi']['Ther']['vec']['Cth_JC'].values
    Rth_DC = para['Swi']['Ther']['vec']['Rth_DC'].values
    Cth_DC = para['Swi']['Ther']['vec']['Cth_DC'].values
    Rth_CA = para['Swi']['Ther']['vec']['Rth_CA'].values
    Cth_CA = para['Swi']['Ther']['vec']['Cth_CA'].values
    Rth_JC_cap = para['Cap']['Ther']['vec']['Rth_JC'].values
    Cth_JC_cap = para['Cap']['Ther']['vec']['Cth_JC'].values
    Rth_CA_cap = para['Cap']['Ther']['vec']['Rth_CA'].values
    Cth_CA_cap = para['Cap']['Ther']['vec']['Cth_CA'].values

    ###################################################################################################################
    # Pre-Processing
    ###################################################################################################################
    Rth_JC = Rth_JC[np.logical_not(pd.isna(Rth_JC))]
    Cth_JC = Cth_JC[np.logical_not(pd.isna(Cth_JC))]
    Rth_DC = Rth_DC[np.logical_not(pd.isna(Rth_DC))]
    Cth_DC = Cth_DC[np.logical_not(pd.isna(Cth_DC))]
    Rth_CA = Rth_CA[np.logical_not(pd.isna(Rth_CA))]
    Cth_CA = Cth_CA[np.logical_not(pd.isna(Cth_CA))]
    Rth_JC_cap = Rth_JC_cap[np.logical_not(pd.isna(Rth_JC_cap))]
    Cth_JC_cap = Cth_JC_cap[np.logical_not(pd.isna(Cth_JC_cap))]
    Rth_CA_cap = Rth_CA_cap[np.logical_not(pd.isna(Rth_CA_cap))]
    Cth_CA_cap = Cth_CA_cap[np.logical_not(pd.isna(Cth_CA_cap))]

    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    if setupPara['Ther']['Heatsink'] == 1: 
        Rth_JA = np.concatenate((Rth_JC, Rth_CA))
        Cth_JA = np.concatenate((Cth_JC, Cth_CA))
        Rth_DA = np.concatenate((Rth_DC, Rth_CA))
        Cth_DA = np.concatenate((Cth_DC, Cth_CA))
        Rth_JA_cap = np.concatenate((Rth_JC_cap, Rth_CA_cap))
        Cth_JA_cap = np.concatenate((Cth_JC_cap, Cth_CA_cap))
    else:
        Rth_JA = Rth_JC
        Cth_JA = Cth_JC
        Rth_DA = Rth_DC
        Cth_DA = Cth_DC
        Rth_JA_cap = Rth_JC_cap
        Cth_JA_cap = Cth_JC_cap

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return [Rth_JA, Cth_JA, Rth_DA, Cth_DA, Rth_JA_cap, Cth_JA_cap]