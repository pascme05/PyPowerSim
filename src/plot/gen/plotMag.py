#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         plotMag
# Date:         12.05.2025
# Author:       Max V. Mueller
# Version:      V.1.0
# Copyright:    Pascal Schirmer
#######################################################################################################################
#######################################################################################################################

#######################################################################################################################
# Function Description
#######################################################################################################################
"""
This function unifies the magnetics/transformer plots for all topologies.
Inputs:     1) time:    time domain results
            2) freq:    frequency domain results
            3) setup:   includes all simulation variables
Outputs:    None
"""

#######################################################################################################################
# Import libs
#######################################################################################################################
# ==============================================================================
# Internal
# ==============================================================================
from src.general.helpFnc import OoM
from src.general.helpFnc import thd

# ==============================================================================
# External
# ==============================================================================
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import matplotlib.gridspec as gridspec
from scipy.fft import fft
import numpy as np
import math
import matplotlib
matplotlib.use('TkAgg')


#######################################################################################################################
# Function
#######################################################################################################################
def plotMag(time, freq, setup, para):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("START: Plotting Magnetics Results")

    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Init
    # ==============================================================================
    timeSw = time['Sw']
    timeAc = time['Ac']
    timeDc = time['Dc']
    freqSw = freq['Sw']
    freqAc = freq['Ac']
    freqDc = freq['Dc']

    # ==============================================================================
    # Parameters
    # ==============================================================================
    fel = setup['Top']['fel']
    fs = setup['Par']['PWM']['fs']
    fsim = setup['Exp']['fsim']
    Q = int(fs / fel)
    R = setup['Top']['R']
    L = setup['Top']['L']
    Mi = setup['Dat']['stat']['Mi']
    Vdc = setup['Dat']['stat']['Vdc']
    phiE = setup['Top']['phiE']
    down = int(setup['Dat']['stat']['cyc']) - 2
    down2 = int(fsim / fs / 200)
    if down2 < 1:
        down2 = 1

    # ==============================================================================
    # Variables
    # ==============================================================================
    t = timeSw['t'].values
    f = fsim * np.linspace(0, 0.5, int(len(t) / 2)) / fel

    ###################################################################################################################
    # Pre-Processing
    ###################################################################################################################
    # ==============================================================================
    # Start and End Plotting
    # ==============================================================================
    # ------------------------------------------
    # Limits
    # ------------------------------------------
    K = int(np.round((t[-1] - t[0]) * fel))
    start = int((len(t) - 1) / K) * (K - 1)
    ende = len(t)

    """
    # ------------------------------------------
    # Change time
    # ------------------------------------------
    t = t[start:ende]
    timeSw = timeSw[:][start:ende]
    for c1 in timeAc:
        timeAc[c1] = timeAc[c1][start:ende]
    for c1 in timeDc:
        timeDc[c1] = timeDc[c1][start:ende]
    """

    t = t[start:ende]
    # ------------------------------------------
    # Transformer currents and voltages
    # ------------------------------------------
    plt.figure()
    txt = "Transformer waveforms for turns ratio of " + "$n_1/n_2$=" + str(para['Tra']['Elec']['con']['n1']/para['Tra']['Elec']['con']['n2'])
    plt.suptitle(txt, size=18)
    plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)

    # Voltages
    plt.subplot(2, 2, 1)
    plt.plot(t[::down2], timeAc['v_1'][::down2], 'b')
    plt.plot(t[::down2], timeAc['v_w1'][::down2], 'b--')
    plt.plot(t[::down2], timeAc['v_2'][::down2], 'r')
    plt.plot(t[::down2], timeAc['v_w2'][::down2], 'r--')
    plt.ylabel("$v(t)$ (V)")
    plt.title('Transformer terminal voltages')
    plt.xlabel('time in (sec)')
    plt.legend(["$v_{1}(t)$", "$v_{w1}(t)$", "$v_{2}(t)$", "$v_{w2}(t)$"], loc='upper right')
    plt.grid('on')

    # Currents
    ax = plt.subplot(2, 2, 2)
    plt.plot(t[::down2], timeAc['i_1'][::down2], 'b')
    plt.plot(t[::down2], timeAc['i_w1'][::down2], 'b--')
    plt.plot(t[::down2], timeAc['i_2'][::down2], 'r')
    plt.plot(t[::down2], timeAc['i_w2'][::down2], 'r--')
    plt.ylabel("$i(t)$ (A)")
    plt.title('Transformer (winding-) currents')
    plt.xlabel('time in (sec)')
    plt.legend(["$i_{1}(t)$", "$i_{w1}(t)$", "$i_{2}(t)$", "$i_{w2}(t)$"], loc='upper right')
    plt.grid('on')

    # Flux density
    ax = plt.subplot(2, 2, 3)
    plt.plot(t[::down2], timeAc['B'][::down2], 'r')
    plt.ylabel("$B(t)$ (T)")
    plt.title('Transformer core flux density')
    plt.xlabel('time in (sec)')
    plt.legend(["$B(t)$"], loc='upper right')
    plt.grid('on')

    ax = plt.subplot(2, 2, 4)
    plt.plot(t[::down2], time['Ther']['tra']['core'][::down2], 'g')
    plt.plot(t[::down2], time['Ther']['tra']['pri'][::down2], 'b')
    plt.plot(t[::down2], time['Ther']['tra']['sec'][::down2], 'r')
    plt.ylabel("$Temperature$ (°C)")
    plt.title('Transformer temperatures')
    plt.xlabel('time in (sec)')
    plt.legend(["Core temperature", "Primary winding temperature", "Secondary winding temperature"], loc='upper right')
    plt.grid('on')

    plt.show()
    ###################################################################################################################
    # MSG OUT
    ###################################################################################################################
    print("END: Plotting Magnetics Results")
