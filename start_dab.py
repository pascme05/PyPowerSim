#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         start_dab
# Date:         04.02.2026
# Author:       Dr. Pascal A. Schirmer
# Version:      V.1.0
# Copyright:    Pascal Schirmer
#######################################################################################################################
#######################################################################################################################

#######################################################################################################################
# Function Description
#######################################################################################################################
"""
This function starts a DAB simulation based on the provided configuration under \config as well as the parameters
set below. The DAB is a dual active bridge converter based on SPS (Single Phase Shift) assumptions.
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
setup['Exp']['Name'] = "dab_test"                                                                                       # name of the simulation (str)
setup['Exp']['Author'] = "Pascal Schirmer"                                                                               # name of the responsible person (str)
setup['Exp']['debug'] = 0                                                                                                # (0): debug mode de-activated, (1): debug mode activated level-1, (2): debug mode activated level-2

# ------------------------------------------
# Operating Mode
# ------------------------------------------
setup['Exp']['output'] = 'I'                                                                                             # (Mi): modulation index controlled, (V): voltage is controlled, (I): current is controlled, (P): active power is controlled, (Q): reactive power is controlled, (Phi): phase shift is controlled (DAB)
setup['Exp']['type'] = 1                                                                                                 # (0): sweep analysis, (1): steady-state analysis, (2): transient analysis, (3): closed loop analysis

# ==============================================================================
# Input Files
# ==============================================================================
# ------------------------------------------
# Mission Profile and Config
# ------------------------------------------
setup['Exp']['conf'] = "default_DCDC"                                                                                    # name of the configuration file

# ------------------------------------------
# Devices
# ------------------------------------------
# Inverter Topologies (B6, B4, B2)
setup['Exp']['Swi'] = "IKQ75N120CS6"                                                                                     # filename of the default switch parameter file
setup['Exp']['Cap'] = "Elco"                                                                                             # filename of the parameter set for the DC link capacitor

# DCDC Topologies (DAB, PSFB)
setup['Exp']['SwiPri'] = "A2F12M12W2"                                                                                    # filename of the primary bridge switch parameter file
setup['Exp']['SwiSec'] = "A2F12M12W2"                                                                                    # filename of the secondary bridge switch parameter file
setup['Exp']['CapPri'] = "none"                                                                                          # filename of the parameter set for the Input capacitor
setup['Exp']['CapSec'] = "none"                                                                                          # filename of the parameter set for the Output capacitor
setup['Exp']['Trafo'] = "trafoDAB"                                                                                       # filename of the parameter set for the transformer

# ==============================================================================
# Plotting and Saving
# ==============================================================================
setup['Exp']['plot'] = 3                                                                                                 # (0): no results are plotted, (1): results are plotted, (2): analytic results are plotted, (3): topology specific plots
setup['Exp']['plotGen'] = 0                                                                                              # (0): no generic plots, (1): loss and thermal models are plotted
setup['Exp']['save'] = 0                                                                                                 # (0): no results are saved, (1): results are saved

#######################################################################################################################
# Calculations
#######################################################################################################################
if __name__ == '__main__':
    main(setup, setupPath)
