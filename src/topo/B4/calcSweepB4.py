#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         calcSweepB4
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
from src.topo.B4.calcSSeqB4 import calcSSeqB4_CB, calcSSeqB4_FF
from src.topo.B4.calcDistB4 import calcDistB4_Ana
from src.general.calcDistNum import calcDistNum
from src.topo.B4.calcTimeB4 import calcTimeB4
from src.general.calcFreq import calcFreq
from src.topo.B4.initB4 import initB4
from src.general.genWaveform import genWave
from src.topo.B4.outB4 import outB4_Sweep

# ==============================================================================
# External
# ==============================================================================
import numpy as np
import math
from tqdm import tqdm


#######################################################################################################################
# Function
#######################################################################################################################
def calcSweepB4(mdl, _, setupTopo, setupData, setupPara, setupExp):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("START: Sweeping B4 bridge")

    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Init
    # ==============================================================================
    v_ref = {}

    # ==============================================================================
    # Parameters
    # ==============================================================================
    fel = setupTopo['fel']
    fsim = setupExp['fsim']
    N = int(fsim / fel)
    K = setupData['stat']['cyc']
    W = setupData['stat']['W']
    Mi = setupData['stat']['Mi']

    # ==============================================================================
    # Variables
    # ==============================================================================
    # ------------------------------------------
    # Init 
    # ------------------------------------------
    [_, _, _, _, _, _, _, distAc, distDc] = initB4(W)

    # ------------------------------------------
    # Inputs 
    # ------------------------------------------
    E = setupTopo['E']
    Vdc = setupData['stat']['Vdc']
    phiE = math.radians(setupTopo['phiE'])
    phiV = math.radians(setupData['stat']['phi'])

    ###################################################################################################################
    # Pre-Processing
    ###################################################################################################################
    # ==============================================================================
    # Generate Reference Waveform
    # ==============================================================================
    # ------------------------------------------
    # Time
    # ------------------------------------------
    t = np.linspace(0, K / fel, K * N + 1)

    # ------------------------------------------
    # Reference
    # ------------------------------------------
    v_ref['A'] = +(Vdc / 2) * Mi * genWave(t, fel, phiV, 0, setupTopo)
    v_ref['B'] = -(Vdc / 2) * Mi * genWave(t, fel, phiV, 0, setupTopo)
    e_ref = E * genWave(t, fel, phiE, 0, setupTopo)

    # ==============================================================================
    # Maximum Modulation Index
    # ==============================================================================
    Mi_max = 1
    if Mi > Mi_max:
        Mi_max = Mi
    M_i = np.linspace(1e-3, Mi_max - 1e-3, W)

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
    if setupPara['PWM']['type'] == "FF":
        [xs, xsh, s, c] = calcSSeqB4_FF(v_ref, t, Mi, setupPara, setupTopo)
    elif setupPara['PWM']['type'] == "CB":
        [xs, xsh, s, c] = calcSSeqB4_CB(v_ref, t, Mi, setupPara, setupTopo)
    else:
        [xs, xsh, s, c] = calcSSeqB4_CB(v_ref, t, Mi, setupPara, setupTopo)

    # ------------------------------------------
    # Time Domain
    # ------------------------------------------
    [timeAc, timeDc] = calcTimeB4(t, s, e_ref, Vdc, Mi, mdl, setupTopo, start, ende)

    # ==============================================================================
    # Sweeping
    # ==============================================================================
    for i in tqdm(range(len(M_i)), desc='Sweep'):
        # ------------------------------------------
        # Switching
        # ------------------------------------------
        if setupPara['PWM']['type'] == "FF":
            [_, _, s, _] = calcSSeqB4_FF(v_ref, t, M_i[i], setupPara, setupTopo)
        elif setupPara['PWM']['type'] == "CB":
            [_, _, s, _] = calcSSeqB4_CB(v_ref, t, M_i[i], setupPara, setupTopo)
        else:
            [_, _, s, _] = calcSSeqB4_CB(v_ref, t, M_i[i], setupPara, setupTopo)

        # ------------------------------------------
        # Time
        # ------------------------------------------
        [tempTimeAc, tempTimeDc] = calcTimeB4(t, s, e_ref, Vdc, M_i[i], mdl, setupTopo, start, ende)

        # ------------------------------------------
        # Distortion
        # ------------------------------------------
        [numDistAc, numDistDc] = calcDistNum(t[start:ende], tempTimeAc['i_a'], tempTimeAc['v_a'],
                                             tempTimeDc['i_dc'], tempTimeDc['v_dc'], 2*Vdc, fel)
        [anaTimeAc, anaTimeDc] = calcDistB4_Ana(M_i[i], Vdc, setupTopo, setupPara)

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
    [freqSw, freqAc, freqDc] = calcFreq(s['A'][start:ende], xs['A'][start:ende], timeAc, timeDc)

    # ==============================================================================
    # Output
    # ==============================================================================
    [time, freq, sweep] = outB4_Sweep(distAc, distDc, timeAc, timeDc, freqSw, freqAc, freqDc, t, v_ref, e_ref, s, c, xs,
                                      xsh, W, M_i, ende, start)

    ###################################################################################################################
    # MSG Out
    ###################################################################################################################
    print("END: Sweeping B4 bridge")

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return [time, freq, sweep]
