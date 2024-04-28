#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         main
# Date:         27.04.2024
# Author:       Dr. Pascal A. Schirmer
# Version:      V.1.0
# Copyright:    Pascal Schirmer
#######################################################################################################################
#######################################################################################################################

#######################################################################################################################
# Function Description
#######################################################################################################################
"""
This function builds the main part of the simulation. It takes the setup files and path variables as input and
executes the program.
Inputs:     1) setup:   includes all simulation variables
            2) path:    includes all path variables
Outputs:    None
"""

#######################################################################################################################
# Import libs
#######################################################################################################################
# ==============================================================================
# Internal
# ==============================================================================
from src.data.loadPara import loadPara
from src.data.loadSetup import loadSetup
from src.topo.B2.calcSweepB2 import calcSweepB2
from src.topo.B4.calcSweepB4 import calcSweepB4
from src.topo.B6.calcSweepB6 import calcSweepB6
from src.topo.B2.calcSteadyB2 import calcSteadyB2
from src.topo.B4.calcSteadyB4 import calcSteadyB4
from src.topo.B6.calcSteadyB6 import calcSteadyB6
from src.topo.B2.calcTransB2 import calcTransB2
from src.topo.B4.calcTransB4 import calcTransB4
from src.topo.B6.calcTransB6 import calcTransB6
from src.plot.plot import plot
from src.plot.plotResults import plotResults
from src.general.genTF import genTF
from src.general.sanityCheck import sanityInput
from src.general.saveResults import saveResults
from src.general.genLoadInput import genLoadInput

# ==============================================================================
# External
# ==============================================================================
import sys


#######################################################################################################################
# Steady-State
#######################################################################################################################
def main(setup, path):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("----------------------------------------------------------------------------------------------------------")
    print("----------------------------------------------------------------------------------------------------------")
    print("Welcome to the Power Electronics Distortion Toolkit!")
    print("Author:     Dr. Pascal Alexander Schirmer")
    print("Copyright:  Pascal Schirmer")
    print("Version:    v.1.0")
    print("Date:       27.04.2024")
    print("----------------------------------------------------------------------------------------------------------")
    print("----------------------------------------------------------------------------------------------------------")

    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Variables
    # ==============================================================================
    time = []
    freq = []
    sweep = []

    ###################################################################################################################
    # Loading
    ###################################################################################################################
    # ==============================================================================
    # MSG IN
    # ==============================================================================
    print("=======================================================================")
    print("START: Loading")
    print("=======================================================================")

    # ==============================================================================
    # Configuration
    # ==============================================================================
    try:
        setup = loadSetup(setup, path)
    except:
        sys.exit('ERROR: Configuration could not be loaded')

    # ==============================================================================
    # Parameter
    # ==============================================================================
    try:
        para = loadPara(setup, path)
    except:
        sys.exit('ERROR: Parameters could not be loaded')

    # ==============================================================================
    # MSG OUT
    # ==============================================================================
    print("=======================================================================")
    print("END: Loading")
    print("=======================================================================")
    print("\n")

    ###################################################################################################################
    # Pre-Processing
    ###################################################################################################################
    # ==============================================================================
    # MSG IN
    # ==============================================================================
    print("=======================================================================")
    print("START: Pre-Processing")
    print("=======================================================================")

    # ==============================================================================
    # Sanity Checks
    # ==============================================================================
    setup = sanityInput(para, setup)

    # ==============================================================================
    # Transfer Functions
    # ==============================================================================
    mdl = genTF(para, setup)

    # ==============================================================================
    # Control Mode
    # ==============================================================================
    setup = genLoadInput(setup)

    # ==============================================================================
    # MSG OUT
    # ==============================================================================
    print("=======================================================================")
    print("END: Pre-Processing")
    print("=======================================================================")
    print("\n")

    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    # ==============================================================================
    # MSG IN
    # ==============================================================================
    print("=======================================================================")
    print("START: Simulation")
    print("=======================================================================")

    # ==============================================================================
    # B2
    # ==============================================================================
    if setup['Top']['sourceType'] == "B2":
        # ------------------------------------------
        # Init Topology
        # ------------------------------------------

        # ------------------------------------------
        # Sweep
        # ------------------------------------------
        if setup['Exp']['type'] == 0:
            [time, freq, sweep] = calcSweepB2(mdl, para, setup)

        # ------------------------------------------
        # Stationary
        # ------------------------------------------
        if setup['Exp']['type'] == 1:
            [time, freq] = calcSteadyB2(mdl, para, setup)

        # ------------------------------------------
        # Transient
        # ------------------------------------------
        if setup['Exp']['type'] == 2:
            [time, freq] = calcTransB2(mdl, para, setup)

    # ==============================================================================
    # B4 
    # ==============================================================================
    elif setup['Top']['sourceType'] == "B4":
        # ------------------------------------------
        # Sweep
        # ------------------------------------------
        if setup['Exp']['type'] == 0:
            [time, freq, sweep] = calcSweepB4(mdl, para, setup)

        # ------------------------------------------
        # Stationary
        # ------------------------------------------
        if setup['Exp']['type'] == 1:
            [time, freq] = calcSteadyB4(mdl, para, setup)

        # ------------------------------------------
        # Transient
        # ------------------------------------------
        if setup['Exp']['type'] == 2:
            [time, freq] = calcTransB4(mdl, para, setup)

    # ==============================================================================
    # B6
    # ==============================================================================
    elif setup['Top']['sourceType'] == "B6":
        # ------------------------------------------
        # Sweep
        # ------------------------------------------
        if setup['Exp']['type'] == 0:
            [time, freq, sweep] = calcSweepB6(mdl, para, setup)

        # ------------------------------------------
        # Stationary
        # ------------------------------------------
        if setup['Exp']['type'] == 1:
            [time, freq] = calcSteadyB6(mdl, para, setup)

        # ------------------------------------------
        # Transient
        # ------------------------------------------
        if setup['Exp']['type'] == 2:
            [time, freq] = calcTransB6(mdl, para, setup)

    # ==============================================================================
    # Default
    # ==============================================================================
    else:
        print("ERROR: Invalid topology")

    # ==============================================================================
    # MSG OUT
    # ==============================================================================
    print("=======================================================================")
    print("END: Simulation")
    print("=======================================================================")
    print("\n")

    ###################################################################################################################
    # Post-Processing
    ###################################################################################################################
    # ==============================================================================
    # MSG IN
    # ==============================================================================
    print("=======================================================================")
    print("START: Post-Processing")
    print("=======================================================================")

    # ==============================================================================
    # Save Results
    # ==============================================================================
    if setup['Exp']['save'] == 1:
        saveResults(time, freq, sweep, setup, path)

    # ==============================================================================
    # Plot Results
    # ==============================================================================
    if setup['Exp']['type'] != 0:
        plotResults(time, setup)

    # ==============================================================================
    # Plot Results
    # ==============================================================================
    if setup['Exp']['plot'] != 0:
        plot(mdl, para, time, freq, sweep, setup)

    # ==============================================================================
    # MSG OUT
    # ==============================================================================
    print("=======================================================================")
    print("END: Post-Processing")
    print("=======================================================================")

    ###################################################################################################################
    # MSG Out
    ###################################################################################################################
