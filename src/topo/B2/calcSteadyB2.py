#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         calcSteadyB2
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
from src.topo.B2.calcSSeqB2 import calcSSeqB2_CB, calcSSeqB2_FF
from src.topo.B2.calcTimeB2 import calcTimeB2
from src.general.calcFreq import calcFreq
from src.elec.calcElecSwi import calcElecSwi
from src.elec.calcLossSwi import calcLossSwi
from src.elec.calcLossCap import calcLossCap
from src.therm.calcTherRC import calcTherRC
from src.topo.B2.initB2 import initB2
from src.general.genWaveform import genWave
from src.therm.initRC import initRC
from src.elec.calcElecCap import calcElecCap
from src.topo.B2.outB2 import outB2_Steady
from src.general.calcAvg import calcAvg

# ==============================================================================
# External
# ==============================================================================
import numpy as np
import math


#######################################################################################################################
# Function
#######################################################################################################################
def calcSteadyB2(mdl, para, setupTopo, setupData, setupPara, setupExp):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("------------------------------------------")
    print("START: Steady-state analysis B2 bridge")
    print("------------------------------------------")

    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Init
    # ==============================================================================
    # ------------------------------------------
    # Variables
    # ------------------------------------------
    out = {}

    # ==============================================================================
    # Parameters
    # ==============================================================================
    # ------------------------------------------
    # General
    # ------------------------------------------
    fel = setupTopo['fel']
    fsim = setupExp['fsim']
    fs = setupPara['PWM']['fs']
    Nsim = int(np.ceil(fsim / fel))
    Npwm = int(np.ceil(fs / fel))
    N = int(fsim / fel)
    K = setupData['stat']['cyc']
    W = setupData['stat']['W']
    Mi = setupData['stat']['Mi']
    start = int(N) * 2
    ende = int(K * N + 1)

    # ------------------------------------------
    # Electrical
    # ------------------------------------------
    E = setupTopo['E']
    Vdc = setupData['stat']['Vdc']
    phiE = math.radians(setupTopo['phiE'])
    phiV = math.radians(setupData['stat']['phi'])

    # ------------------------------------------
    # Thermal
    # ------------------------------------------
    Tc = setupData['stat']['Tc']

    # ------------------------------------------
    # Solver
    # ------------------------------------------
    err = np.inf
    T_sw = np.ones((3, 1)) * setupData['stat']['Tj']
    T_ca = setupData['stat']['Tj']
    T_old = setupData['stat']['Tj']
    iter = 0

    # ==============================================================================
    # Variables
    # ==============================================================================
    [_, timeElec, timeLoss, timeTher, _, _, _, _, _] = initB2(W)
    
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
    v_ref = (Vdc/2) * Mi * genWave(t, fel, phiV, 0, setupTopo)
    e_ref = E * genWave(t, fel, phiE, 0, setupTopo) 

    # ==============================================================================
    # Thermal ROM
    # ==============================================================================
    # ------------------------------------------
    # Load
    # ------------------------------------------
    [Rth_JA, Cth_JA, Rth_DA, Cth_DA, Rth_CA, Cth_CA, Rth_JA_cap, Cth_JA_cap] = initRC(para, setupPara)

    # ------------------------------------------
    # Variables
    # ------------------------------------------
    Tinit_T = np.zeros((len(Rth_JA), 2))
    Tinit_D = np.zeros((len(Rth_DA), 2))
    Tinit_K = np.zeros((len(Rth_CA), 2))
    Tinit_C = np.zeros(np.size(Rth_JA_cap))

    ###################################################################################################################
    # Calculation (Stationary)
    ###################################################################################################################
    # ==============================================================================
    # Switching Function
    # ==============================================================================
    if setupPara['PWM']['type'] == "FF":
        [xs, xsh, s, c] = calcSSeqB2_FF(v_ref, t, Mi, setupPara, setupTopo)
    elif setupPara['PWM']['type'] == "CB":
        [xs, xsh, s, c] = calcSSeqB2_CB(v_ref, t, Mi, setupPara, setupTopo)
    else:
        [xs, xsh, s, c] = calcSSeqB2_CB(v_ref, t, Mi, setupPara, setupTopo)

    # ==============================================================================
    # Time Domain
    # ==============================================================================
    [timeAc, timeDc] = calcTimeB2(t, s, e_ref, Vdc, Mi, mdl, setupTopo, start, ende)

    # ==============================================================================
    # Msg
    # ==============================================================================
    print("START: Iterating")
    print("------------------------------------------")

    # ==============================================================================
    # Iterating
    # ==============================================================================
    while err > setupExp['tol']:
        # ------------------------------------------
        # Electrical
        # ------------------------------------------
        # Switches
        timeElec['sw']['S1'] = calcElecSwi(Vdc, timeAc['i_a'], (s[start:ende] == +1), T_sw[0], 'HS', para, setupPara)
        timeElec['sw']['S2'] = calcElecSwi(Vdc, timeAc['i_a'], (s[start:ende] == -1), T_sw[1], 'LS', para, setupPara)

        # Capacitor
        timeDc['v_dc'] = calcElecCap(t, timeDc['i_c'], T_ca, para, setupPara, setupTopo)
        timeElec['cap']['C1']['i_c'] = timeDc['i_c']
        timeElec['cap']['C1']['v_c'] = timeDc['v_dc']

        # ------------------------------------------
        # Losses
        # ------------------------------------------
        # Switches
        timeLoss['sw']['S1'] = calcLossSwi(s[start:ende] * (+1), timeElec['sw']['S1']['i_T'],
                                           timeElec['sw']['S1']['i_D'], timeElec['sw']['S1']['v_T'],
                                           timeElec['sw']['S1']['v_D'], T_sw[0], para, setupPara, setupExp)
        timeLoss['sw']['S2'] = calcLossSwi(s[start:ende] * (-1), timeElec['sw']['S2']['i_T'],
                                           timeElec['sw']['S2']['i_D'], timeElec['sw']['S2']['v_T'],
                                           timeElec['sw']['S2']['v_D'], T_sw[1], para, setupPara, setupExp)

        # Capacitor
        timeLoss['cap']['C1'] = calcLossCap(t, timeDc['i_c'], T_ca, para, setupPara, setupTopo)

        # ------------------------------------------
        # Init Thermal
        # ------------------------------------------
        if iter == 0:
            # Switches
            Tinit_T[:, 0] = np.mean(timeLoss['sw']['S1']['p_T']) * Rth_JA
            Tinit_D[:, 0] = np.mean(timeLoss['sw']['S1']['p_D']) * Rth_DA
            Tinit_T[:, 1] = np.mean(timeLoss['sw']['S2']['p_T']) * Rth_JA
            Tinit_D[:, 1] = np.mean(timeLoss['sw']['S2']['p_D']) * Rth_DA
            if setupPara['Ther']['Coupling'] == 1:
                Tinit_K[:, 0] = np.mean(timeLoss['sw']['S1']['p_L']) * Rth_CA
                Tinit_K[:, 1] = np.mean(timeLoss['sw']['S2']['p_L']) * Rth_CA

            # Capacitor
            Tinit_C = np.mean(timeLoss['cap']['C1']['p_L']) * Rth_JA_cap

        # ------------------------------------------
        # Thermal
        # ------------------------------------------
        # Switches
        [timeTher['sw']['T1'], Tinit_T[:, 0]] = calcTherRC(Tinit_T[:, 0], Tc, timeLoss['sw']['S1']['p_T'], t[start:ende], Rth_JA, Cth_JA)
        [timeTher['sw']['T2'], Tinit_T[:, 1]] = calcTherRC(Tinit_T[:, 1], Tc, timeLoss['sw']['S2']['p_T'], t[start:ende], Rth_JA, Cth_JA)
        [timeTher['sw']['D1'], Tinit_D[:, 0]] = calcTherRC(Tinit_D[:, 0], Tc, timeLoss['sw']['S1']['p_D'], t[start:ende], Rth_DA, Cth_DA)
        [timeTher['sw']['D2'], Tinit_D[:, 1]] = calcTherRC(Tinit_D[:, 1], Tc, timeLoss['sw']['S2']['p_D'], t[start:ende], Rth_DA, Cth_DA)

        if setupPara['Ther']['Coupling'] == 1:
            [timeTher['sw']['C1'], Tinit_K[:, 0]] = calcTherRC(Tinit_K[:, 0], Tc, timeLoss['sw']['S1']['p_L'], t[start:ende], Rth_CA, Cth_CA)
            [timeTher['sw']['C2'], Tinit_K[:, 1]] = calcTherRC(Tinit_K[:, 1], Tc, timeLoss['sw']['S1']['p_L'], t[start:ende], Rth_CA, Cth_CA)
            timeTher['sw']['T1'] = timeTher['sw']['T1'][:] + timeTher['sw']['C1'] - Tc
            timeTher['sw']['D1'] = timeTher['sw']['D1'][:] + timeTher['sw']['C1'] - Tc
            timeTher['sw']['T2'] = timeTher['sw']['T2'][:] + timeTher['sw']['C2'] - Tc
            timeTher['sw']['D2'] = timeTher['sw']['D2'][:] + timeTher['sw']['C2'] - Tc
        else:
            timeTher['sw']['C1'] = Tc * np.ones(np.size(s[start:ende]))
            timeTher['sw']['C2'] = Tc * np.ones(np.size(s[start:ende]))

        # Capacitor
        [timeTher['cap']['C1'], Tinit_C] = calcTherRC(Tinit_C, Tc, timeLoss['cap']['C1']['p_L'], t[start:ende], Rth_JA_cap, Cth_JA_cap)

        # ------------------------------------------
        # Error
        # ------------------------------------------
        err = np.mean(np.abs(timeTher['sw']['T1'] - T_old)) / np.mean(T_old)

        # ------------------------------------------
        # Update Para
        # ------------------------------------------
        if setupExp['therFeed'] == 1:
            # Switches
            T_sw[0] = np.mean(timeTher['sw']['T1'])
            T_sw[1] = np.mean(timeTher['sw']['T2'])

            # Capacitor
            T_ca = np.mean(timeTher['cap']['C1'])

        # Iteration
        iter = iter + 1

        # Previous Temperature
        T_old = timeTher['sw']['T1']

        # ------------------------------------------
        # Msg
        # ------------------------------------------
        if iter < int(setupExp['int']):
            print("ITER: %d) Stationary temperature T_swi=%.2f C (T_cap=%.2f C) and P_swi=%.2f W (Pv_cap=%.2f W) with error: %.2f %%" % (
                  iter, T_sw[0], T_ca, np.mean(timeLoss['sw']['S1']['p_L']), np.mean(timeLoss['cap']['C1']['p_L']), err * 100))
        else:
            print("ITER: %d) Maximum iteration reached" % iter)
            break

    # ==============================================================================
    # Msg
    # ==============================================================================
    print("------------------------------------------")
    print("INFO: Converged after %d iterations with T_swi=%.2f C (T_cap=%.2f C) and error: %.2f %%" % (iter, T_sw[0], T_ca, err * 100))

    ###################################################################################################################
    # Post-Processing
    ###################################################################################################################
    # ==============================================================================
    # Averaging
    # ==============================================================================
    # ------------------------------------------
    # Adapt
    # ------------------------------------------
    out['elec'] = timeElec
    out['loss'] = timeLoss

    # ------------------------------------------
    # Average
    # ------------------------------------------
    out = calcAvg(out, setupExp, setupTopo, Nsim, Npwm)

    # ------------------------------------------
    # Removing
    # ------------------------------------------
    for c0 in out:
        for c1 in out[c0]:
            for c2 in out[c0][c1]:
                out[c0][c1][c2] = out[c0][c1][c2].iloc[0:int(ende - start)]
                out[c0][c1][c2] = out[c0][c1][c2].reset_index(drop=True)

    # ------------------------------------------
    # Out
    # ------------------------------------------
    timeElec = out['elec']
    timeLoss = out['loss']

    # ==============================================================================
    # Update Thermal
    # ==============================================================================
    # ------------------------------------------
    # Switches
    # ------------------------------------------
    # Normal
    [timeTher['sw']['T1'], _] = calcTherRC(Tinit_T[:, 0], Tc, timeLoss['sw']['S1']['p_T'], t[start:ende], Rth_JA, Cth_JA)
    [timeTher['sw']['T2'], _] = calcTherRC(Tinit_T[:, 1], Tc, timeLoss['sw']['S2']['p_T'], t[start:ende], Rth_JA, Cth_JA)
    [timeTher['sw']['D1'], _] = calcTherRC(Tinit_D[:, 0], Tc, timeLoss['sw']['S1']['p_D'], t[start:ende], Rth_DA, Cth_DA)
    [timeTher['sw']['D2'], _] = calcTherRC(Tinit_D[:, 1], Tc, timeLoss['sw']['S2']['p_D'], t[start:ende], Rth_DA, Cth_DA)

    # Coupling
    if setupPara['Ther']['Coupling'] == 1:
        [timeTher['sw']['C1'], _] = calcTherRC(Tinit_K[:, 0], Tc, timeLoss['sw']['S1']['p_L'], t[start:ende], Rth_CA, Cth_CA)
        [timeTher['sw']['C2'], _] = calcTherRC(Tinit_K[:, 1], Tc, timeLoss['sw']['S1']['p_L'], t[start:ende], Rth_CA, Cth_CA)
        timeTher['sw']['T1'] = timeTher['sw']['T1'][:] + timeTher['sw']['C1'] - Tc
        timeTher['sw']['D1'] = timeTher['sw']['D1'][:] + timeTher['sw']['C1'] - Tc
        timeTher['sw']['T2'] = timeTher['sw']['T2'][:] + timeTher['sw']['C2'] - Tc
        timeTher['sw']['D2'] = timeTher['sw']['D2'][:] + timeTher['sw']['C2'] - Tc
    else:
        timeTher['sw']['C1'] = Tc * np.ones(np.size(s[start:ende]))
        timeTher['sw']['C2'] = Tc * np.ones(np.size(s[start:ende]))
    # ------------------------------------------
    # Capacitor
    # ------------------------------------------
    [timeTher['cap']['C1'], _] = calcTherRC(Tinit_C, Tc, timeLoss['cap']['C1']['p_L'], t[start:ende], Rth_JA_cap, Cth_JA_cap)

    # ==============================================================================
    # Frequency domain
    # ==============================================================================
    [freqSw, freqAc, freqDc] = calcFreq(s[start:ende], xs[start:ende], timeAc, timeDc)
    
    # ==============================================================================
    # Output
    # ==============================================================================
    [time, freq] = outB2_Steady(timeElec, timeLoss, timeTher, timeAc, timeDc, freqSw, freqAc, freqDc, t, v_ref, e_ref,
                                s, c, xs, xsh, W, ende, start)
    
    ###################################################################################################################
    # MSG Out
    ###################################################################################################################
    print("------------------------------------------")
    print("END: Steady-state analysis B2 bridge")
    print("------------------------------------------")

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return [time, freq]
