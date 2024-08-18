#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         calcSweep
# Date:         04.05.2024
# Author:       Dr. Pascal A. Schirmer
# Version:      V.1.0
# Copyright:    Pascal Schirmer
#######################################################################################################################
#######################################################################################################################

#######################################################################################################################
# Function Description
#######################################################################################################################
"""
This function calculates the results for a parameter sweep for any given topology class.
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
from src.general.calcDistNum import calcDistNum

# ==============================================================================
# External
# ==============================================================================
import numpy as np
import math
from tqdm import tqdm


#######################################################################################################################
# Function
#######################################################################################################################
def calcSweep(top, mdl, _, setup):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("------------------------------------------")
    print("START: Sweeping class", top.name)
    print("------------------------------------------")

    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Parameters
    # ==============================================================================
    fel = setup['Top']['fel']
    fsim = setup['Exp']['fsim']
    N = int(fsim / fel)
    K = int(setup['Dat']['stat']['cyc'])
    W = int(setup['Dat']['stat']['W'])
    Mi = setup['Dat']['stat']['Mi']
    Mi_max = top.Mi_max
    E = setup['Top']['E']
    Vdc = setup['Dat']['stat']['Vdc']
    phiE = math.radians(setup['Top']['phiE'])
    phiV = math.radians(setup['Dat']['stat']['phi'])
    L = setup['Top']['L']
    R = setup['Top']['R']
    Z = np.sqrt(R ** 2 + (2 * np.pi * fel * L) ** 2)

    # ==============================================================================
    # Variables
    # ==============================================================================
    [_, _, _, _, _, _, _, distAc, distDc] = top.initOut()

    ###################################################################################################################
    # Pre-Processing
    ###################################################################################################################
    # ==============================================================================
    # Generate Reference Waveform
    # ==============================================================================
    t_ref = np.linspace(0, K / fel, K * N + 1)
    [v_ref, e_ref, _] = top.calcRef(E, phiE, phiV, [], setup)

    # ==============================================================================
    # Maximum Modulation Index
    # ==============================================================================
    if Mi > Mi_max:
        Mi = Mi_max
    M_i = np.linspace(1e-3, Mi - 1e-3, W)

    # ==============================================================================
    # Start and End
    # ==============================================================================
    start = int(N) * 2
    ende = int(K * N + 1)

    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    # ==============================================================================
    # Stationary
    # ==============================================================================
    # ------------------------------------------
    # Switching Function
    # ------------------------------------------
    [xs, xsh, s, c, x, xN0] = top.calcPWM(v_ref, t_ref, Mi, setup)

    # ------------------------------------------
    # Time Domain
    # ------------------------------------------
    [timeAc, timeDc, _] = top.calcTime(s, e_ref, t_ref, Mi, mdl, start, ende, [], 1, setup)

    # ==============================================================================
    # Sweeping
    # ==============================================================================
    for i in tqdm(range(len(M_i)), desc='Sweep'):
        # ------------------------------------------
        # Switching
        # ------------------------------------------
        [_, _, s_i, _, _, _] = top.calcPWM(v_ref, t_ref, M_i[i], setup)

        # ------------------------------------------
        # Time
        # ------------------------------------------
        [tempAc, tempDc, _] = top.calcTime(s_i, e_ref, t_ref, M_i[i], mdl, start, ende, [], 1, setup)

        # ------------------------------------------
        # Distortion
        # ------------------------------------------
        [numDistAc, numDistDc] = calcDistNum(t_ref[start:ende], tempAc['i_a'], tempAc['v_a'], tempDc['i_dc'], tempDc['v_dc'], Vdc, fel)
        [anaTimeAc, anaTimeDc] = top.calcDist(tempAc['i_a'], tempAc['v_a'], M_i[i], L, Z, setup)

        # ------------------------------------------
        # Output
        # ------------------------------------------
        for c1 in numDistAc:
            distAc['num'][c1][i] = numDistAc[c1]
            distAc['ana'][c1][i] = anaTimeAc[c1]
        for c1 in numDistDc:
            distDc['num'][c1][i] = numDistDc[c1]
            distDc['ana'][c1][i] = anaTimeDc[c1]

    ###################################################################################################################
    # Post-Processing
    ###################################################################################################################
    # ==============================================================================
    # Frequency domain
    # ==============================================================================
    [freqSw, freqAc, freqDc] = calcFreq(s['A'][start:ende], xs['A'][start:ende], timeAc['i_a'], timeAc['v_a'],
                                        timeAc['v_a0'], timeDc['i_dc'], timeDc['v_dc'])

    # ==============================================================================
    # Output
    # ==============================================================================
    [time, freq, sweep] = top.out([], [], [], timeAc, timeDc, freqSw, freqAc, freqDc, distAc, distDc, t_ref, v_ref,
                                  e_ref, s, c, xs, xsh, x, xN0, M_i, start, ende, 1)

    ###################################################################################################################
    # MSG Out
    ###################################################################################################################
    print("------------------------------------------")
    print("END: Sweeping class", top.name)
    print("------------------------------------------")

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return [time, freq, sweep]
