#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         calcTherRC
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


#######################################################################################################################
# Function
#######################################################################################################################
def calcTherRC(Tinit, Tc, Pv, t, Rth, Cth):
    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Parameters
    # ==============================================================================
    N = len(Pv)
    K = len(Rth)
    tau = Rth*Cth

    # ==============================================================================
    # Variables
    # ==============================================================================
    dt = np.diff(t)
    dt = np.insert(dt, len(dt), dt)
    T = Tinit*np.ones((N, K))

    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    for i in range(1, N):
        T[i, :] = (2*tau[:]-dt[i])/(2*tau[:]+dt[i])*T[i-1, :] + (Rth[:]*dt[i])/(2*tau[:]+dt[i])*(Pv[i]+Pv[i-1])

    ###################################################################################################################
    # Post-Processing
    ###################################################################################################################
    Tj = np.sum(T, axis=1) + Tc
    Tend = T[-1, :]

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return [Tj, Tend]
