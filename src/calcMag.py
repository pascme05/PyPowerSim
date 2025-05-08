#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         calcMag
# Date:         08.05.2025
# Author:       Max V. Mueller
# Version:      V.1.0
# Copyright:    Pascal Schirmer
#######################################################################################################################
#######################################################################################################################

#######################################################################################################################
# Function Description
#######################################################################################################################
"""
This function calculates the results for a converter including a transformer as a magnetic element, e.g. a full bridge
connected to a transformer, as in a Dual-Active-Bridge (DAB) or Phase-Shifted-Full-Bridge (PSFB) converter
Inputs:     1) top:     topology class
            2) mdl:     all models and transfer functions of the architecture
            3) para:    all parameters used in the simulation
            4) setup:   includes all simulation variables
Outputs:    1) time:    results in the time domain
            2) freq:    results in the frequency domain
            3) dist:    results in the distortion domain
"""

#######################################################################################################################
# Import libs
#######################################################################################################################
# ==============================================================================
# Internal
# ==============================================================================
from src.general.calcFreq import calcFreq
from src.elec.calcElecSwi import calcElecSwi
from src.elec.calcLossSwi import calcLossSwi
from src.elec.calcLossCap import calcLossCap
from src.therm.calcTherRC import calcTherRC
from src.therm.initRC import initRC
from src.elec.calcElecCap import calcElecCap
from src.general.calcAvg import calcAvg

# ==============================================================================
# External
# ==============================================================================
import numpy as np
import math


#######################################################################################################################
# Function
#######################################################################################################################
def calcMag(top, mdl, para, setup):
    print("Hello world!")
