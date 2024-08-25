#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         journalSteadySecVC
# Date:         08.05.2024
# Author:       Dr. Pascal A. Schirmer
# Version:      V.1.0
# Copyright:    Pascal Schirmer
#######################################################################################################################
#######################################################################################################################

#######################################################################################################################
# Function Description
#######################################################################################################################
"""
This function start a simulation based on the provided configuration under \config as well as the parameters defined in
the start script.
"""

#######################################################################################################################
# Import external libs
#######################################################################################################################
# ==============================================================================
# Internal
# ==============================================================================
from src.main import main
from src.general.helpFnc import initSetup, initPath

# ==============================================================================
# External
# ==============================================================================
import warnings

#######################################################################################################################
# Format
#######################################################################################################################
warnings.filterwarnings("ignore")

#######################################################################################################################
# Paths
#######################################################################################################################
setupPath = initPath('PyPowerSim')

#######################################################################################################################
# Init
#######################################################################################################################
setup = initSetup()

#######################################################################################################################
# Configuration
#######################################################################################################################
# ==============================================================================
# Experiment
# ==============================================================================
# ------------------------------------------
# General
# ------------------------------------------
setup['Exp']['Name'] = "journalSteadySecVC"                                                                              # name of the simulation (str)
setup['Exp']['Author'] = "Pascal Schirmer"                                                                               # name of the responsible person (str)
setup['Exp']['debug'] = 0                                                                                                # (0): debug mode de-activated, (1): debug mode activated level-1, (2): debug mode activated level-2

# ------------------------------------------
# Operating Mode
# ------------------------------------------
setup['Exp']['output'] = 'Mi'                                                                                            # (Mi): modulation index controlled, (V): voltage is controlled, (I): current is controlled, (P): active power is controlled, (Q): reactive power is controlled
setup['Exp']['type'] = 1                                                                                                 # (0): sweep analysis, (1): steady-state analysis, (2): transient analysis, (3): closed loop analysis

# ==============================================================================
# Input Files
# ==============================================================================
# ------------------------------------------
# Mission Profile and Config
# ------------------------------------------
setup['Exp']['conf'] = "journalSteadySecVC"

# ------------------------------------------
# Devices
# ------------------------------------------
setup['Exp']['Swi'] = "2MBI300XBE120"                                                                                    # filename of the parameter set for the switching devices
setup['Exp']['Cap'] = "Elco"                                                                                             # filename of the parameter set for the DC link capacitor

# ==============================================================================
# Plotting and Saving
# ==============================================================================
setup['Exp']['plot'] = 2                                                                                                 # (0): no results are plotted, (1): results are plotted, (2): analytic results are plotted
setup['Exp']['plotGen'] = 1                                                                                              # (0): no generic plots, (1): loss and thermal models are plotted
setup['Exp']['save'] = 0                                                                                                 # (0): no results are saved, (1): results are saved

#######################################################################################################################
# Calculations
#######################################################################################################################
if __name__ == '__main__':
    main(setup, setupPath)
