#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         calcSteadyB6
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
    [_, timeElec, timeLoss, timeTher, _, _, _, _, _] = initB6(W)

    # ------------------------------------------
    # Inputs 
    # ------------------------------------------
    # Electrical
    E = setupTopo['E']
    Vdc = setupData['stat']['Vdc']
    phiE = math.radians(setupTopo['phiE'])
    phiV = math.radians(setupData['stat']['phi'])

    # Thermal
    Tj = setupData['stat']['Tj']
    Tcap = setupData['stat']['Tj']
    Tc = setupData['stat']['Tc']

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
    # Start and End
    # ==============================================================================
    start = int(N) * 2
    ende = int(K * N + 1)

    # ==============================================================================
    # Thermal ROM
    # ==============================================================================
    [Rth_JA, Cth_JA, Rth_DA, Cth_DA, Rth_CA, Cth_CA, Rth_JA_cap, Cth_JA_cap] = initRC(para, setupPara)

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
    # Open-Loop Losses
    # ==============================================================================
    if setupExp['loop'] == "OL":
        # ------------------------------------------
        # Electrical
        # ------------------------------------------
        # Switches
        for i in range(0, len(id2)):
            timeElec['sw'][id2[i]] = calcElecSwi(Vdc, timeAc[id4[i]], (s[id3[i]][start:ende] == (-1) ** i), Tj, id5[i],
                                                 para, setupPara)

        # Capacitor
        timeDc['v_dc'] = calcElecCap(t, timeDc['i_c'], Tcap, para, setupPara, setupTopo)
        timeElec['cap']['C1']['i_c'] = timeDc['i_c']
        timeElec['cap']['C1']['v_c'] = timeDc['v_dc']

        # ------------------------------------------
        # Losses
        # ------------------------------------------
        # Switches
        for i in range(0, len(id2)):
            timeLoss['sw'][id2[i]] = calcLossSwi(s[id3[i]][start:ende] * (-1) ** i, timeElec['sw'][id2[i]]['i_T'],
                                                 timeElec['sw'][id2[i]]['i_D'], timeElec['sw'][id2[i]]['v_T'],
                                                 timeElec['sw'][id2[i]]['v_D'], Tj, para, setupPara, setupExp)

        # Capacitor
        timeLoss['cap']['C1'] = calcLossCap(t, timeDc['i_c'], Tcap, para, setupPara, setupTopo)

    # ==============================================================================
    # Closed-Loop Losses
    # ==============================================================================
    if setupExp['loop'] == "CL":
        # ------------------------------------------
        # Init
        # ------------------------------------------
        err = np.inf
        T_old = np.ones((7, 1)) * Tj
        T_new = np.ones((7, 1)) * Tj
        iter = 0

        # ------------------------------------------
        # Iterating
        # ------------------------------------------
        # Msg
        print("START: Iterating")
        print("------------------------------------------")

        # Start
        while err > setupExp['tol']:
            # Electrical
            for i in range(0, len(id2)):
                timeElec['sw'][id2[i]] = calcElecSwi(Vdc, timeAc[id4[i]], (s[id3[i]][start:ende] == (-1) ** i),
                                                     T_old[i], id5[i], para, setupPara)

            timeDc['v_dc'] = calcElecCap(t, timeDc['i_c'], T_old[6], para, setupPara, setupTopo)
            timeElec['cap']['C1']['i_c'] = timeDc['i_c']
            timeElec['cap']['C1']['v_c'] = timeDc['v_dc']

            # Losses
            for i in range(0, len(id2)):
                timeLoss['sw'][id2[i]] = calcLossSwi(s[id3[i]][start:ende] * (-1) ** i, timeElec['sw'][id2[i]]['i_T'],
                                                     timeElec['sw'][id2[i]]['i_D'], timeElec['sw'][id2[i]]['v_T'],
                                                     timeElec['sw'][id2[i]]['v_D'], T_old[i], para, setupPara, setupExp)
            timeLoss['cap']['C1'] = calcLossCap(t, timeDc['i_c'], T_old[6], para, setupPara, setupTopo)

            # Thermal
            for i in range(0, len(id2)):
                if setupPara['Ther']['Heatsink'] == 1 & setupPara['Ther']['Coupling'] == 1:
                    T_new[i] = np.mean(timeLoss['sw'][id2[i]]['p_T']) * np.sum(Rth_JA) + np.mean(
                        timeLoss['sw'][id2[i]]['p_L']) * np.sum(Rth_CA) + Tc
                else:
                    T_new[i] = np.mean(timeLoss['sw'][id2[i]]['p_T']) * np.sum(Rth_JA) + Tc
            T_new[6] = np.mean(timeLoss['cap']['C1']['p_L']) * np.sum(Rth_JA_cap) + Tc

            # Error
            err = np.sum(abs(T_old - T_new)) / np.sum(T_old)

            # Update
            T_old = T_new[:]
            iter = iter + 1

            # Msg
            print("ITER: %d) Stationary temperature T_swi=%.2f (T_cap=%.2f) and P_swi=%.2f (Pv_cap=%.2f) with error: %.3f" % (
                  iter, T_old[0], T_old[6], np.mean(timeLoss['sw']['S1']['p_L']), np.mean(timeLoss['cap']['C1']['p_L']), err * 100))

        # ------------------------------------------
        # Converged
        # ------------------------------------------
        # Msg
        print("------------------------------------------")
        print("INFO: Converged after %d iterations with T_swi=%.2f (T_cap=%.2f) and error: %.3f" % (iter, T_old[0], T_old[6], err * 100))

    # ==============================================================================
    # Transient Thermal
    # ==============================================================================
    # ------------------------------------------
    # Msg
    # ------------------------------------------
    print("INFO: Calculating transient temperatures")

    # ------------------------------------------
    # Thermal
    # ------------------------------------------
    # Switches
    for i in range(0, len(id2)):
        [timeTher['sw'][id6[i]], _] = calcTherRC(0, Tc, timeLoss['sw'][id2[i]]['p_T'], t[start:ende], Rth_JA, Cth_JA)
        [timeTher['sw'][id7[i]], _] = calcTherRC(0, Tc, timeLoss['sw'][id2[i]]['p_D'], t[start:ende], Rth_DA, Cth_DA)

    # Capacitor
    [timeTher['cap']['C1'], _] = calcTherRC(0, Tc, timeLoss['cap']['C1']['p_L'], t[start:ende], Rth_JA_cap, Cth_JA_cap)

    # Coupling
    for i in range(0, len(id2)):
        if setupPara['Ther']['Heatsink'] == 1 & setupPara['Ther']['Coupling'] == 1:
            [timeTher['sw'][id8[i]], _] = calcTherRC(0, Tc, timeLoss['sw'][id2[i]]['p_L'], t[start:ende], Rth_CA, Cth_CA)
            timeTher['sw'][id6[i]] = timeTher['sw'][id6[i]][:] + timeTher['sw'][id8[i]][:] - Tc
            timeTher['sw'][id7[i]] = timeTher['sw'][id7[i]][:] + timeTher['sw'][id8[i]][:] - Tc
        else:
            timeTher['sw'][id8[i]] = Tc * np.ones(np.size(s['A'][start:ende]))

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
    [time, freq] = outB6_Steady(timeElec, timeLoss, timeTher, timeAc, timeDc, freqSw, freqAc, freqDc, t, v_ref, e_ref,
                                s, c, xs, xsh, x, n0, W, ende, start)

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
