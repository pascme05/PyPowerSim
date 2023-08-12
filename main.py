#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         main
# Date:         01.14.2023
# Author:       Dr. Pascal A. Schirmer
# Version:      V.0.1
# Copyright:    Pascal Schirmer
#######################################################################################################################
#######################################################################################################################

#######################################################################################################################
# Import libs
#######################################################################################################################
# ==============================================================================
# Internal
# ==============================================================================
from src.data.loadPara import loadPara
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


#######################################################################################################################
# Steady-State
#######################################################################################################################
def main(setupExp, setupData, setupTopo, setupPara, setupPath):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("----------------------------------------------------------------------------------------------------------")
    print("----------------------------------------------------------------------------------------------------------")
    print("Welcome to the Power Electronics Distortion Toolkit!")
    print("Author:     Dr. Pascal Alexander Schirmer")
    print("Copyright:  Pascal Schirmer")
    print("Version:    v.0.1")
    print("Date:       01.04.2023")
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
    para = []

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
    # Parameter
    # ==============================================================================
    try:
        para = loadPara(setupTopo, setupPath, setupPara)
    except:
        print("ERROR: Parameters could not be loaded")

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
    [setupExp, setupData, setupTopo, setupPara] = sanityInput(setupExp, setupData, setupTopo, setupPara)
    
    # ==============================================================================
    # Transfer Functions
    # ==============================================================================
    mdl = genTF(para, setupTopo)

    # ==============================================================================
    # Control Mode
    # ==============================================================================
    setupData = genLoadInput(setupExp, setupTopo, setupData)

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
    if setupTopo['sourceType'] == "B2":
        # ------------------------------------------
        # Sweep
        # ------------------------------------------
        if setupExp['type'] == 0:
            [time, freq, sweep] = calcSweepB2(mdl, para, setupTopo, setupData, setupPara, setupExp)

        # ------------------------------------------
        # Stationary
        # ------------------------------------------
        if setupExp['type'] == 1:
            [time, freq] = calcSteadyB2(mdl, para, setupTopo, setupData, setupPara, setupExp)

        # ------------------------------------------
        # Transient
        # ------------------------------------------
        if setupExp['type'] == 2:
            [time, freq] = calcTransB2(mdl, para, setupTopo, setupData, setupPara, setupExp)
        

    # ==============================================================================
    # B4 
    # ==============================================================================
    elif setupTopo['sourceType'] == "B4":
        # ------------------------------------------
        # Sweep
        # ------------------------------------------
        if setupExp['type'] == 0:
            [time, freq, sweep] = calcSweepB4(mdl, para, setupTopo, setupData, setupPara, setupExp)

        # ------------------------------------------
        # Stationary
        # ------------------------------------------
        if setupExp['type'] == 1:
            [time, freq] = calcSteadyB4(mdl, para, setupTopo, setupData, setupPara, setupExp)

        # ------------------------------------------
        # Transient
        # ------------------------------------------
        if setupExp['type'] == 2:
            [time, freq] = calcTransB4(mdl, para, setupTopo, setupData, setupPara, setupExp)

    # ==============================================================================
    # B6
    # ==============================================================================
    elif setupTopo['sourceType'] == "B6":
        # ------------------------------------------
        # Sweep
        # ------------------------------------------
        if setupExp['type'] == 0:
            [time, freq, sweep] = calcSweepB6(mdl, para, setupTopo, setupData, setupPara, setupExp)

        # ------------------------------------------
        # Stationary
        # ------------------------------------------
        if setupExp['type'] == 1:
            [time, freq] = calcSteadyB6(mdl, para, setupTopo, setupData, setupPara, setupExp)

        # ------------------------------------------
        # Transient
        # ------------------------------------------
        if setupExp['type'] == 2:
            [time, freq] = calcTransB6(mdl, para, setupTopo, setupData, setupPara, setupExp)

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
    if setupExp['save'] == 1:
        saveResults(time, freq, sweep, setupExp, setupData, setupPara, setupTopo, setupPath)

    # ==============================================================================
    # Plot Results
    # ==============================================================================
    if setupExp['type'] != 0:
        plotResults(time, setupTopo)

    # ==============================================================================
    # Plot Results
    # ==============================================================================
    plot(time, freq, sweep, setupPara, setupData, setupTopo, setupExp)
    
    # ==============================================================================
    # MSG OUT
    # ==============================================================================
    print("=======================================================================")
    print("END: Post-Processing")
    print("=======================================================================")

    ###################################################################################################################
    # MSG Out
    ###################################################################################################################