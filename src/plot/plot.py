#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         plot
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
from src.plot.spe.plotSweep_B2 import plotSweep_B2
from src.plot.spe.plotStat_B2 import plotStat_B2
from src.plot.spe.plotTrans_B2 import plotTrans_B2
from src.plot.spe.plotSweep_B4 import plotSweep_B4
from src.plot.spe.plotStat_B4 import plotStat_B4
from src.plot.spe.plotTrans_B4 import plotTrans_B4
from src.plot.spe.plotSweep_B6 import plotSweep_B6
from src.plot.spe.plotStat_B6 import plotStat_B6
from src.plot.spe.plotTrans_B6 import plotTrans_B6
from src.plot.gen.plotGen import plotGenTF, plotGenLoss, plotGenTher
from src.plot.gen.plotSweep import plotSweep
from src.plot.gen.plotStat import plotStat
from src.plot.gen.plotTrans import plotTrans
from src.plot.gen.plotClose import plotClose
from src.plot.gen.plotMag import plotMag


# ==============================================================================
# External
# ==============================================================================
import matplotlib.style as mplstyle
mplstyle.use('fast')


#######################################################################################################################
# Function
#######################################################################################################################
def plot(mdl, para, time, freq, sweep, setup):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("------------------------------------------")
    print("START: Plotting")
    print("------------------------------------------")

    ###################################################################################################################
    # General
    ###################################################################################################################
    if setup['Exp']['plotGen'] == 1:
        # ------------------------------------------
        # Msg
        # ------------------------------------------
        print("INFO: Plotting general information")

        # ------------------------------------------
        # Bode Plots
        # ------------------------------------------
        try:
            plotGenTF(mdl, setup['Par'], setup['Top'])
        except:
            print("WARN: Bode plots could not been created")

        # ------------------------------------------
        # Loss Models
        # ------------------------------------------
        try:
            plotGenLoss(para, setup['Par'], setup['Dat'])
        except:
            print("WARN: Loss plots could not been created")

        # ------------------------------------------
        # Thermal Models
        # ------------------------------------------
        try:
            plotGenTher(para)
        except:
            print("WARN: Thermal plots could not been created")



    ###################################################################################################################
    # Topology
    ###################################################################################################################
    # ==============================================================================
    # Generic
    # ==============================================================================
    if setup['Exp']['plot'] != 3:
        # ------------------------------------------
        # Sweeping
        # ------------------------------------------
        if setup['Exp']['type'] == 0:
            plotSweep(time, freq, sweep, setup)

        # ------------------------------------------
        # Steady State
        # ------------------------------------------
        elif setup['Exp']['type'] == 1:
            plotStat(time, freq, setup)

        # ------------------------------------------
        # Transient
        # ------------------------------------------
        elif setup['Exp']['type'] == 2:
            plotTrans(time, freq, setup)

        # ------------------------------------------
        # Closed Loop
        # ------------------------------------------
        elif setup['Exp']['type'] == 3:
            plotClose(time, freq, setup)

    # ==============================================================================
    # Magnetics
    # ==============================================================================
    if setup['Top']['LD_tra'] != 'NT' and setup['Top']['sourceType'] != "B6":
        plotMag(time, freq, setup)

    # ==============================================================================
    # Specific
    # ==============================================================================
    else:
        # ------------------------------------------
        # B2
        # ------------------------------------------
        if setup['Top']['sourceType'] == "B2":
            # Sweep
            if setup['Exp']['type'] == 0:
                plotSweep_B2(time, freq, sweep, setup['Par'], setup['Dat'], setup['Top'], setup['Exp'])

            # Stationary
            if setup['Exp']['type'] == 1:
                plotStat_B2(time, freq, setup['Par'], setup['Dat'], setup['Top'], setup['Exp'])

            # Transient
            if setup['Exp']['type'] == 2:
                plotTrans_B2(time, freq, setup['Par'], setup['Dat'], setup['Top'], setup['Exp'])

        # ------------------------------------------
        # B4
        # ------------------------------------------
        if setup['Top']['sourceType'] == "B4":
            # Sweep
            if setup['Exp']['type'] == 0:
                plotSweep_B4(time, freq, sweep, setup['Par'], setup['Dat'], setup['Top'], setup['Exp'])

            # Stationary
            if setup['Exp']['type'] == 1:
                plotStat_B4(time, freq, setup['Par'], setup['Dat'], setup['Top'], setup['Exp'])

            # Transient
            if setup['Exp']['type'] == 2:
                plotTrans_B4(time, freq, setup['Par'], setup['Dat'], setup['Top'], setup['Exp'])

        # ------------------------------------------
        # B6
        # ------------------------------------------
        if setup['Top']['sourceType'] == "B6":
            # Sweep
            if setup['Exp']['type'] == 0:
                plotSweep_B6(time, freq, sweep, setup['Par'], setup['Dat'], setup['Top'], setup['Exp'])

            # Stationary
            if setup['Exp']['type'] == 1:
                plotStat_B6(time, freq, setup['Par'], setup['Dat'], setup['Top'], setup['Exp'])

            # Transient
            if setup['Exp']['type'] == 2:
                plotTrans_B6(time, freq, setup['Par'], setup['Dat'], setup['Top'], setup['Exp'])

    ###################################################################################################################
    # MSG Out
    ###################################################################################################################
    print("------------------------------------------")
    print("END: Plotting waveforms")
    print("------------------------------------------")

    ###################################################################################################################
    # Return
    ###################################################################################################################
