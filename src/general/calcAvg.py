#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         calcAvg
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
This function calculates the average of power losses assuming fundamental or switching frequency.
Inputs:     1) data:    input data array
            2) setup:   all setup variables
            3) Nsim:    number of simulation samples
            4) Npwm:    number of samples per PWM period
Outputs:    1) data:    output data array
"""

#######################################################################################################################
# Import libs
#######################################################################################################################
# ==============================================================================
# Internal
# ==============================================================================
import pandas as pd

# ==============================================================================
# External
# ==============================================================================
from src.topo.B2.initB2 import initB2_Data
from src.topo.B4.initB4 import initB4_Data
from src.topo.B6.initB6 import initB6_Data


#######################################################################################################################
# Function
#######################################################################################################################
def calcAvg(data, setup, Nsim, Npwm):
    ###################################################################################################################
    # Init
    ###################################################################################################################
    if setup['Top']['sourceType'] == 'B2':
        out = initB2_Data()
    elif setup['Top']['sourceType'] == 'B4':
        out = initB4_Data()
    elif setup['Top']['sourceType'] == 'B6':
        out = initB6_Data()
    else:
        out = initB6_Data()

    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    # ==============================================================================
    # Electrical
    # ==============================================================================
    if setup['Exp']['freqAvg'] == 'fel':
        for c1 in data['elec']:
            for c2 in data['elec'][c1]:
                for c3 in data['elec'][c1][c2]:
                    temp = pd.concat([data['elec'][c1][c2][c3][0:Nsim+1], data['elec'][c1][c2][c3], data['elec'][c1][c2][c3][-1-Nsim:-1]]).rolling(window=Nsim, center=True, closed='both', min_periods=None).mean()
                    out['elec'][c1][c2][c3] = temp[Nsim:-Nsim]
    elif setup['Exp']['freqAvg'] == 'fs':
        for c1 in data['elec']:
            for c2 in data['elec'][c1]:
                for c3 in data['elec'][c1][c2]:
                    temp = pd.concat([data['elec'][c1][c2][c3][0:Nsim+1], data['elec'][c1][c2][c3], data['elec'][c1][c2][c3][-1-Nsim:-1]]).rolling(window=int(Nsim/Npwm), center=True, closed='both', min_periods=None).mean()
                    out['elec'][c1][c2][c3] = temp[Nsim:-Nsim]

    # ==============================================================================
    # Losses
    # ==============================================================================
    if setup['Exp']['freqAvg'] == 'fel':
        for c1 in data['loss']:
            for c2 in data['loss'][c1]:
                for c3 in data['loss'][c1][c2]:
                    temp = pd.concat([data['loss'][c1][c2][c3][0:Nsim+1], data['loss'][c1][c2][c3], data['loss'][c1][c2][c3][-1-Nsim:-1]]).rolling(window=Nsim, center=True, closed='both', min_periods=None).mean()
                    out['loss'][c1][c2][c3] = temp[Nsim:-Nsim]
    elif setup['Exp']['freqAvg'] == 'fs':
        for c1 in data['loss']:
            for c2 in data['loss'][c1]:
                for c3 in data['loss'][c1][c2]:
                    temp = pd.concat([data['loss'][c1][c2][c3][0:Nsim+1], data['loss'][c1][c2][c3], data['loss'][c1][c2][c3][-1-Nsim:-1]]).rolling(window=int(Nsim/Npwm), center=True, closed='both', min_periods=None).mean()
                    out['loss'][c1][c2][c3] = temp[Nsim:-Nsim]

    ###################################################################################################################
    # Post
    ###################################################################################################################
    if setup['Exp']['freqAvg'] == 'fel' or setup['Exp']['freqAvg'] == 'fs':
        data = out
        
    ###################################################################################################################
    # Return
    ###################################################################################################################
    return data
