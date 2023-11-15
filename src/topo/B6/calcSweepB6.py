#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         calcSweepB6
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
from src.topo.B6.calcSSeqB6 import calcSSeqB6_CB, calcSSeqB6_FF, calcSSeqB6_SV, calcSSeqB6_OPP
from src.topo.B6.calcDistB6 import calcDistB6_Ana
from src.topo.B6.calcTimeB6 import calcTimeB6
from src.general.calcFreq import calcFreq
from src.general.calcDistNum import calcDistNum
from src.topo.B6.initB6 import initB6
from src.pwm.genWaveform import genWave
from src.topo.B6.outB6 import outB6_Sweep

# ==============================================================================
# External
# ==============================================================================
import numpy as np
import math
from tqdm import tqdm


#######################################################################################################################
# Function
#######################################################################################################################
def calcSweepB6(mdl, _, setupTopo, setupData, setupPara, setupExp):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("START: Sweeping B6 bridge")

    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Init
    # ==============================================================================
    # ------------------------------------------
    # Variables
    # ------------------------------------------
    v_ref = {}
    e_ref = {}

    # ------------------------------------------
    # IDs
    # ------------------------------------------
    id1 = ['A', 'B', 'C']

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
    [_, _, _, _, _, _, _, distAc, distDc] = initB6(W)

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
    for i in range(0, len(id1)):
        v_ref[id1[i]] = (Vdc / 2) * Mi * genWave(t, fel, phiV, -i * 2 / 3 * np.pi, setupTopo)
        e_ref[id1[i]] = E * genWave(t, fel, phiE, -i * 2 / 3 * np.pi, setupTopo)

    # ==============================================================================
    # Maximum Modulation Index
    # ==============================================================================
    if setupPara['PWM']['zero'] == "SPWM":
        Mi_max = 1.0
    else:
        Mi_max = 2 / np.sqrt(3)
    if Mi > Mi_max:
        Mi_max = Mi
    M_i = np.linspace(0, Mi_max - 1e-3, W)

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
        [xs, xsh, s, c, x, n0] = calcSSeqB6_FF(v_ref, t, Mi, setupPara, setupTopo)
    elif setupPara['PWM']['type'] == "CB":
        [xs, xsh, s, c, x, n0] = calcSSeqB6_CB(v_ref, t, Mi, setupPara, setupTopo)
    elif setupPara['PWM']['type'] == "SV":
        [xs, xsh, s, c, x, n0] = calcSSeqB6_SV(v_ref, t, Mi, setupPara, setupTopo)
    elif setupPara['PWM']['type'] == "OPP":
        [xs, xsh, s, c, x, n0] = calcSSeqB6_OPP(v_ref, t, Mi, setupPara, setupTopo)
    else:
        [xs, xsh, s, c, x, n0] = calcSSeqB6_CB(v_ref, t, Mi, setupPara, setupTopo)

    # ------------------------------------------
    # Time Domain
    # ------------------------------------------
    [timeAc, timeDc] = calcTimeB6(t, s, e_ref, Vdc, Mi, mdl, setupTopo, start, ende)

    # ==============================================================================
    # Sweeping
    # ==============================================================================
    for i in tqdm(range(len(M_i)), desc='Sweep'):
        # ------------------------------------------
        # Switching
        # ------------------------------------------
        if setupPara['PWM']['type'] == "FF":
            [_, _, s_i, _, _, _] = calcSSeqB6_FF(v_ref, t, M_i[i], setupPara, setupTopo)
        elif setupPara['PWM']['type'] == "CB":
            [_, _, s_i, _, _, _] = calcSSeqB6_CB(v_ref, t, M_i[i], setupPara, setupTopo)
        elif setupPara['PWM']['type'] == "SV":
            [_, _, s_i, _, _, _] = calcSSeqB6_SV(v_ref, t, M_i[i], setupPara, setupTopo)
        elif setupPara['PWM']['type'] == "OPP":
            [_, _, s_i, _, _, _] = calcSSeqB6_OPP(v_ref, t, M_i[i], setupPara, setupTopo)
        else:
            [_, _, s_i, _, _, _] = calcSSeqB6_CB(v_ref, t, M_i[i], setupPara, setupTopo)

        # ------------------------------------------
        # Time
        # ------------------------------------------
        [tempTimeAc, tempTimeDc] = calcTimeB6(t, s_i, e_ref, Vdc, M_i[i], mdl, setupTopo, start, ende)

        # ------------------------------------------
        # Distortion
        # ------------------------------------------
        [numDistAc, numDistDc] = calcDistNum(t[start:ende], tempTimeAc['i_a'], tempTimeAc['v_a'], tempTimeDc['i_dc'],
                                             tempTimeDc['v_dc'], Vdc, fel)
        [anaTimeAc, anaTimeDc] = calcDistB6_Ana(t[start:ende], tempTimeAc['i_a'], tempTimeAc['v_a'],
                                                numDistAc['I_a_v1_eff'], M_i[i], Vdc, setupTopo, setupPara)

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
    [time, freq, sweep] = outB6_Sweep(distAc, distDc, timeAc, timeDc, freqSw, freqAc, freqDc, t, v_ref, e_ref, s, c, xs,
                                      xsh, x, n0, W, M_i, ende, start)

    ###################################################################################################################
    # MSG Out
    ###################################################################################################################
    print("END: Sweeping B6 bridge")

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return [time, freq, sweep]
