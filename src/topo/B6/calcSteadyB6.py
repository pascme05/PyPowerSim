#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         calcSteadyB6
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
from src.topo.B6.calcSSeqB6 import calcSSeqB6_CB, calcSSeqB6_FF, calcSSeqB6_SV
from src.topo.B6.calcTimeB6 import calcTimeB6
from src.general.calcFreq import calcFreq
from src.elec.calcElecSwi import calcElecSwi
from src.elec.calcLossSwi import calcLossSwi
from src.elec.calcLossCap import calcLossCap
from src.therm.calcTherRC import calcTherRC
from src.topo.B6.initB6 import initB6
from src.general.genWaveform import genWave
from src.therm.initRC import initRC
from src.elec.calcElecCap import calcElecCap
from src.topo.B6.outB6 import outB6_Steady
from src.general.calcAvg import calcAvg

# ==============================================================================
# External
# ==============================================================================
import numpy as np
import math


#######################################################################################################################
# Function
#######################################################################################################################
def calcSteadyB6(mdl, para, setupTopo, setupData, setupPara, setupExp):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("------------------------------------------")
    print("START: Steady-state analysis B6 bridge")
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
    v_ref = {}
    e_ref = {}
    out = {}

    # ------------------------------------------
    # IDs
    # ------------------------------------------
    id1 = ['A', 'B', 'C']
    id2 = ['S1', 'S2', 'S3', 'S4', 'S5', 'S6']
    id3 = ['A', 'A', 'B', 'B', 'C', 'C']
    id4 = ['i_a', 'i_a', 'i_b', 'i_b', 'i_c', 'i_c']
    id5 = ['HS', 'LS', 'HS', 'LS', 'HS', 'LS']
    id6 = ['T1', 'T2', 'T3', 'T4', 'T5', 'T6']
    id7 = ['D1', 'D2', 'D3', 'D4', 'D5', 'D6']
    id8 = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6']

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
    T_sw = np.ones((7, 1)) * setupData['stat']['Tj']
    T_ca = setupData['stat']['Tj']
    T_old = setupData['stat']['Tj']
    iter = 0

    # ==============================================================================
    # Variables
    # ==============================================================================
    [_, timeElec, timeLoss, timeTher, _, _, _, _, _] = initB6(W)

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
    # Thermal ROM
    # ==============================================================================
    # ------------------------------------------
    # Load
    # ------------------------------------------
    [Rth_JA, Cth_JA, Rth_DA, Cth_DA, Rth_CA, Cth_CA, Rth_JA_cap, Cth_JA_cap] = initRC(para, setupPara)

    # ------------------------------------------
    # Variables
    # ------------------------------------------
    Tinit_T = np.zeros((len(Rth_JA), len(id2)))
    Tinit_D = np.zeros((len(Rth_DA), len(id2)))
    Tinit_K = np.zeros((len(Rth_CA), len(id2)))
    Tinit_C = np.zeros(np.size(Rth_JA_cap))

    ###################################################################################################################
    # Calculation (Stationary)
    ###################################################################################################################
    # ==============================================================================
    # Switching Function
    # ==============================================================================
    if setupPara['PWM']['type'] == "FF":
        [xs, xsh, s, c, x, n0] = calcSSeqB6_FF(v_ref, t, Mi, setupPara, setupTopo)
    elif setupPara['PWM']['type'] == "CB":
        [xs, xsh, s, c, x, n0] = calcSSeqB6_CB(v_ref, t, Mi, setupPara, setupTopo)
    elif setupPara['PWM']['type'] == "SV":
        [xs, xsh, s, c, x, n0] = calcSSeqB6_SV(v_ref, t, Mi, setupPara, setupTopo)
    else:
        [xs, xsh, s, c, x, n0] = calcSSeqB6_CB(v_ref, t, Mi, setupPara, setupTopo)

    # ==============================================================================
    # Time Domain
    # ==============================================================================
    [timeAc, timeDc] = calcTimeB6(t, s, e_ref, Vdc, Mi, mdl, setupTopo, start, ende)

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
        for i in range(0, len(id2)):
            timeElec['sw'][id2[i]] = calcElecSwi(Vdc, timeAc[id4[i]], (s[id3[i]][start:ende] == (-1) ** i), T_sw[i],
                                                 id5[i], para, setupPara)

        # Capacitor
        timeDc['v_dc'] = calcElecCap(t, timeDc['i_c'], T_ca, para, setupPara, setupTopo)
        timeElec['cap']['C1']['i_c'] = timeDc['i_c']
        timeElec['cap']['C1']['v_c'] = timeDc['v_dc']

        # ------------------------------------------
        # Losses
        # ------------------------------------------
        # Switches
        for i in range(0, len(id2)):
            timeLoss['sw'][id2[i]] = calcLossSwi(s[id3[i]][start:ende] * (-1) ** i, timeElec['sw'][id2[i]]['i_T'],
                                                 timeElec['sw'][id2[i]]['i_D'], timeElec['sw'][id2[i]]['v_T'],
                                                 timeElec['sw'][id2[i]]['v_D'], T_sw[i], para, setupPara, setupExp)

        # Capacitor
        timeLoss['cap']['C1'] = calcLossCap(t, timeDc['i_c'], T_ca, para, setupPara, setupTopo)

        # ------------------------------------------
        # Init Thermal
        # ------------------------------------------
        if iter == 0:
            # Switches
            for i in range(0, len(id2)):
                Tinit_T[:, i] = np.mean(timeLoss['sw'][id2[i]]['p_T']) * Rth_JA
                Tinit_D[:, i] = np.mean(timeLoss['sw'][id2[i]]['p_D']) * Rth_DA
                if setupPara['Ther']['Coupling'] == 1:
                    Tinit_K[:, i] = np.mean(timeLoss['sw'][id2[i]]['p_L']) * Rth_CA

            # Capacitor
            Tinit_C = np.mean(timeLoss['cap']['C1']['p_L']) * Rth_JA_cap

        # ------------------------------------------
        # Thermal
        # ------------------------------------------
        # Switches
        for i in range(0, len(id2)):
            [timeTher['sw'][id6[i]], Tinit_T[:, i]] = calcTherRC(Tinit_T[:, i], Tc, timeLoss['sw'][id2[i]]['p_T'], t[start:ende], Rth_JA, Cth_JA)
            [timeTher['sw'][id7[i]], Tinit_D[:, i]] = calcTherRC(Tinit_D[:, i], Tc, timeLoss['sw'][id2[i]]['p_D'], t[start:ende], Rth_DA, Cth_DA)
            if setupPara['Ther']['Coupling'] == 1:
                [timeTher['sw'][id8[i]], Tinit_K[:, i]] = calcTherRC(Tinit_K[:, i], Tc, timeLoss['sw'][id2[i]]['p_L'], t[start:ende], Rth_CA, Cth_CA)
                timeTher['sw'][id6[i]] = timeTher['sw'][id6[i]][:] + timeTher['sw'][id8[i]][:] - Tc
                timeTher['sw'][id7[i]] = timeTher['sw'][id7[i]][:] + timeTher['sw'][id8[i]][:] - Tc
            else:
                timeTher['sw'][id8[i]] = Tc * np.ones(np.size(s['A'][start:ende]))

        # Capacitor
        [timeTher['cap']['C1'], Tinit_C] = calcTherRC(Tinit_C, Tc, timeLoss['cap']['C1']['p_L'], t[start:ende], Rth_JA_cap, Cth_JA_cap)

        # ------------------------------------------
        # Error
        # ------------------------------------------
        err = np.mean(np.abs(timeTher['sw']['T1'] - T_old)) / np.mean(T_old)

        # ------------------------------------------
        # Update Para
        # ------------------------------------------
        if setupExp['loop'] == "CL":
            # Switches
            for i in range(0, len(id2)):
                T_sw[i] = np.mean(timeTher['sw'][id6[i]])

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
                  iter, T_sw[0], T_ca, np.mean(timeLoss['sw']['S1']['p_L']), np.mean(timeLoss['cap']['C1']['p_L']), err*100))
        else:
            print("ITER: %d) Maximum iteration reached" % iter)
            break

    # ==============================================================================
    # Msg
    # ==============================================================================
    print("------------------------------------------")
    print("INFO: Converged after %d iterations with T_swi=%.2f C (T_cap=%.2f C) and error: %.2f %%" % (iter, T_sw[0], T_ca, err*100))

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
                out[c0][c1][c2] = out[c0][c1][c2].iloc[0:int(ende-start)]
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
    for i in range(0, len(id2)):
        [timeTher['sw'][id6[i]], _] = calcTherRC(Tinit_T[:, i], Tc, timeLoss['sw'][id2[i]]['p_T'], t[start:ende], Rth_JA, Cth_JA)
        [timeTher['sw'][id7[i]], _] = calcTherRC(Tinit_D[:, i], Tc, timeLoss['sw'][id2[i]]['p_D'], t[start:ende], Rth_DA, Cth_DA)
        if setupPara['Ther']['Coupling'] == 1:
            [timeTher['sw'][id8[i]], _] = calcTherRC(Tinit_K[:, i], Tc, timeLoss['sw'][id2[i]]['p_L'], t[start:ende], Rth_CA, Cth_CA)
            timeTher['sw'][id6[i]] = timeTher['sw'][id6[i]][:] + timeTher['sw'][id8[i]][:] - Tc
            timeTher['sw'][id7[i]] = timeTher['sw'][id7[i]][:] + timeTher['sw'][id8[i]][:] - Tc
        else:
            timeTher['sw'][id8[i]] = Tc * np.ones(np.size(s['A'][start:ende]))

    # ------------------------------------------
    # Capacitor
    # ------------------------------------------
    [timeTher['cap']['C1'], _] = calcTherRC(Tinit_C, Tc, timeLoss['cap']['C1']['p_L'], t[start:ende], Rth_JA_cap, Cth_JA_cap)

    # ==============================================================================
    # Frequency domain
    # ==============================================================================
    [freqSw, freqAc, freqDc] = calcFreq(s['A'][start:ende], xs['A'][start:ende], timeAc, timeDc)

    # ==============================================================================
    # Output
    # ==============================================================================
    [time, freq] = outB6_Steady(timeElec, timeLoss, timeTher, timeAc, timeDc, freqSw, freqAc, freqDc, t, v_ref,
                                e_ref, s, c, xs, xsh, x, n0, W, ende, start)

    ###################################################################################################################
    # MSG Out
    ###################################################################################################################
    print("------------------------------------------")
    print("END: Steady-state analysis B6 bridge")
    print("------------------------------------------")

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return [time, freq]
