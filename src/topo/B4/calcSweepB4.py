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
from src.topo.B4.calcDistB4 import calcDistB4_Num, calcDistB4_Ana 
from src.topo.B4.calcTimeB4 import calcTimeB4
from src.topo.B4.calcFreqB4 import calcFreqB4
from src.topo.B4.initB4 import initB4
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
def calcSweepB4(mdl, para, setupTopo, setupData, setupPara, setupExp):
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
    [timeSw, _, _, _, freqSw, freqDc, freqAc, distAc, distDc] = initB4(W)

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
    v_ref['A'] = +(Vdc/2) * Mi * genWave(t, fel, phiV, 0, setupTopo)
    v_ref['B'] = -(Vdc/2) * Mi * genWave(t, fel, phiV, 0, setupTopo)
    e_ref = E * genWave(t, fel, phiE, 0, setupTopo) 

    # ==============================================================================
    # Maximum Modulation Index
    # ==============================================================================
    Mi_max = 1
    if Mi > Mi_max:
        Mi_max = Mi
    M_i = np.linspace(1e-3, Mi_max-1e-3, W)
    
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
        [xs, xsh, s, c] = calcSSeqB4_FF(v_ref, t, Mi, setupPara, setupTopo)
    elif setupPara['PWM']['type'] == "CB":
        [xs, xsh, s, c] = calcSSeqB4_CB(v_ref, t, Mi, setupPara, setupTopo)

    # ------------------------------------------
    # Time Domain
    # ------------------------------------------
    [timeAc, timeDc] = calcTimeB4(t, s, e_ref, Vdc, Mi, mdl, setupTopo, start, ende)

    # ==============================================================================
    # Sweeping
    # ==============================================================================
    for i in trange(len(M_i), desc='Sweep'):
        # ------------------------------------------
        # Switching
        # ------------------------------------------
        if setupPara['PWM']['type'] == "FF":
            [_, _, s, _] = calcSSeqB4_FF(v_ref, t, M_i[i], setupPara, setupTopo)
        elif setupPara['PWM']['type'] == "CB":
            [_, _, s, _] = calcSSeqB4_CB(v_ref, t, M_i[i], setupPara, setupTopo)
        
        # ------------------------------------------
        # Time
        # ------------------------------------------
        [tempTimeAc, tempTimeDc] = calcTimeB4(t, s, e_ref, Vdc, M_i[i], mdl, setupTopo, start, ende)
        
        # ------------------------------------------
        # Distortion
        # ------------------------------------------
        [numDistAc, numDistDc] = calcDistB4_Num(t[start:ende], tempTimeAc['i_a'], tempTimeAc['v_ab'], tempTimeDc['i_dc'], tempTimeDc['v_dc'], Vdc, setupTopo)
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
    [freqSw, freqAc, freqDc] = calcFreqB4(s['A'][start:ende], xs['A'][start:ende], timeAc, timeDc)
    
    # ==============================================================================
    # Output
    # ==============================================================================
    timeSw['t'] = t[0:(ende-start)]
    timeSw['v_a_ref'] = v_ref['A'][start:ende]
    timeSw['v_b_ref'] = v_ref['B'][start:ende]
    timeSw['e'] = e_ref[start:ende]
    timeSw['sA'] = s['A'][start:ende]
    timeSw['sB'] = s['B'][start:ende]
    timeSw['cA'] = c['A'][start:ende]
    timeSw['cB'] = c['B'][start:ende]
    timeSw['xAs'] = xs['A'][start:ende]
    timeSw['xBs'] = xs['B'][start:ende]
    timeSw['xAsh'] = xsh['A'][start:ende]
    timeSw['xBsh'] = xsh['B'][start:ende]
    
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
    print("END: Sweeping B4 bridge")

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return [time, freq, sweep]