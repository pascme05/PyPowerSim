#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         calcSweepB6
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
from src.topo.B6.calcSSeqB6 import calcSSeqB6_CB, calcSSeqB6_FF, calcSSeqB6_SV
from src.topo.B6.calcDistB6 import calcDistB6_Num, calcDistB6_Ana 
from src.topo.B6.calcTimeB6 import calcTimeB6
from src.topo.B6.calcFreqB6 import calcFreqB6
from src.topo.B6.initB6 import initB6
from src.general.genWaveform import genWave

# ==============================================================================
# External
# ==============================================================================
import numpy as np
import math
from tqdm.auto import trange

#######################################################################################################################
# Function
#######################################################################################################################
def calcSweepB6(mdl, para, setupTopo, setupData, setupPara, setupExp):
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
    v_ref = {}
    e_ref = {}

    # ==============================================================================
    # Parameters
    # ==============================================================================
    fel = setupTopo['fel']
    fsim = setupExp['fsim']
    N = int(fsim/fel)
    K = setupData['stat']['cyc']*2
    W = setupData['stat']['W'] 
    Mi = setupData['stat']['Mi']

    # ==============================================================================
    # Variables
    # ==============================================================================
    # ------------------------------------------
    # Init 
    # ------------------------------------------
    [timeSw, _, _, _, freqSw, freqDc, freqAc, distAc, distDc] = initB6(W)

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
    t = np.linspace(0, K/fel, K*N+1)
    v_ref['A'] = (Vdc/2) * Mi * genWave(t, fel, phiV, 0, setupTopo)
    v_ref['B'] = (Vdc/2) * Mi * genWave(t, fel, phiV, -2/3*np.pi, setupTopo)
    v_ref['C'] = (Vdc/2) * Mi * genWave(t, fel, phiV, -4/3*np.pi, setupTopo)
    e_ref['A'] = E * genWave(t, fel, phiE, 0, setupTopo) 
    e_ref['B'] = E * genWave(t, fel, phiE, -2/3*np.pi, setupTopo) 
    e_ref['C'] = E * genWave(t, fel, phiE, -4/3*np.pi, setupTopo) 

    # ==============================================================================
    # Maximum Modulation Index
    # ==============================================================================
    if setupPara['PWM']['zero'] == "SPWM":
        Mi_max = 1.0
    else:
        Mi_max = 2/np.sqrt(3)
    if Mi > Mi_max:
        Mi_max = Mi
    M_i = np.linspace(0, Mi_max-1e-3, W)
    
    # ==============================================================================
    # Start and End
    # ==============================================================================
    start = int(N * (K/2))
    ende = int(K*N + 1)
    
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

    # ------------------------------------------
    # Time Domain
    # ------------------------------------------
    [timeAc, timeDc] = calcTimeB6(t, s, e_ref, Vdc, Mi, mdl, setupTopo, start, ende)

    # ==============================================================================
    # Sweeping
    # ==============================================================================
    for i in trange(len(M_i), desc='Sweep'):
        # ------------------------------------------
        # Switching
        # ------------------------------------------
        if setupPara['PWM']['type'] == "FF":
            [_, _, s, _, _, _] = calcSSeqB6_FF(v_ref, t, M_i[i], setupPara, setupTopo)
        elif setupPara['PWM']['type'] == "CB":
            [_, _, s, _, _, _] = calcSSeqB6_CB(v_ref, t, M_i[i], setupPara, setupTopo)
        elif setupPara['PWM']['type'] == "SV":
            [_, _, s, _, _, _] = calcSSeqB6_SV(v_ref, t, M_i[i], setupPara, setupTopo)
        
        # ------------------------------------------
        # Time
        # ------------------------------------------
        [tempTimeAc, tempTimeDc] = calcTimeB6(t, s, e_ref, Vdc, M_i[i], mdl, setupTopo, start, ende)
        
        # ------------------------------------------
        # Distortion
        # ------------------------------------------
        [numDistAc, numDistDc] = calcDistB6_Num(t[start:ende], tempTimeAc['i_a'], tempTimeAc['v_a'], tempTimeDc['i_dc'], tempTimeDc['v_dc'], Vdc, setupTopo)
        [anaTimeAc, anaTimeDc] = calcDistB6_Ana(t[start:ende], tempTimeAc['i_a'], tempTimeAc['v_a'], numDistAc['I_a_v1_eff'], M_i[i], Vdc, setupTopo, setupPara)
        
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
    [freqSw, freqAc, freqDc] = calcFreqB6(s['A'][start:ende], xs['A'][start:ende], timeAc, timeDc)
    
    # ==============================================================================
    # Output
    # ==============================================================================
    timeSw['t'] = t[0:(ende-start)]
    timeSw['v_a_ref'] = v_ref['A'][start:ende]
    timeSw['v_b_ref'] = v_ref['B'][start:ende]
    timeSw['v_c_ref'] = v_ref['C'][start:ende]
    timeSw['e_a'] = e_ref['A'][start:ende]
    timeSw['e_b'] = e_ref['B'][start:ende]
    timeSw['e_c'] = e_ref['C'][start:ende]
    timeSw['sA'] = s['A'][start:ende]
    timeSw['sB'] = s['B'][start:ende]
    timeSw['sC'] = s['C'][start:ende]
    timeSw['c'] = c[start:ende]
    timeSw['xAs'] = xs['A'][start:ende]
    timeSw['xBs'] = xs['B'][start:ende]
    timeSw['xCs'] = xs['C'][start:ende]
    timeSw['xAsh'] = xsh['A'][start:ende]
    timeSw['xBsh'] = xsh['B'][start:ende]
    timeSw['xCsh'] = xsh['C'][start:ende]
    timeSw['xA'] = x['A'][start:ende]
    timeSw['xB'] = x['B'][start:ende]
    timeSw['xC'] = x['C'][start:ende]
    timeSw['n0'] = n0[start:ende]
    
    # ==============================================================================
    # Combine
    # ==============================================================================
    time = {}
    time['Sw'] = timeSw
    time['Ac'] = timeAc
    time['Dc'] = timeDc
    freq = {}
    freq['Sw'] = freqSw
    freq['Ac'] = freqAc
    freq['Dc'] = freqDc
    sweep = {}
    sweep['Ac'] = distAc
    sweep['Dc'] = distDc
    sweep['Mi'] = M_i
    
    ###################################################################################################################
    # MSG Out
    ###################################################################################################################
    print("END: Sweeping B6 bridge")

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return [time, freq, sweep]