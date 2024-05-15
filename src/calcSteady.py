#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         calcSteady
# Date:         07.05.2024
# Author:       Dr. Pascal A. Schirmer
# Version:      V.1.0
# Copyright:    Pascal Schirmer
#######################################################################################################################
#######################################################################################################################

#######################################################################################################################
# Function Description
#######################################################################################################################
"""
This function calculates the results for a steady state solution for any given topology class.
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
from src.elec.calcElecSwi import calcElecSwi
from src.elec.calcLossSwi import calcLossSwi
from src.elec.calcLossCap import calcLossCap
from src.therm.calcTherRC import calcTherRC
from src.therm.initRC import initRC
from src.elec.calcElecCap import calcElecCap
from src.general.calcAvg import calcAvg

# ==============================================================================
# External
# ==============================================================================
import numpy as np
import math


#######################################################################################################################
# Function
#######################################################################################################################
def calcSteady(top, mdl, para, setup):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("------------------------------------------")
    print("START: Steady-State solution class", top.name)
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
    fel = setup['Top']['fel']
    fsim = setup['Exp']['fsim']
    N = int(fsim / fel)
    K = int(setup['Dat']['stat']['cyc'])
    Mi = setup['Dat']['stat']['Mi']
    start = int(N) * 2
    ende = int(K * N + 1)

    # ------------------------------------------
    # Electrical
    # ------------------------------------------
    E = setup['Top']['E']
    Vdc = setup['Dat']['stat']['Vdc']
    phiE = math.radians(setup['Top']['phiE'])
    phiV = math.radians(setup['Dat']['stat']['phi'])

    # ------------------------------------------
    # Thermal
    # ------------------------------------------
    Tc = setup['Dat']['stat']['Tc']

    # ------------------------------------------
    # Solver
    # ------------------------------------------
    err = np.inf
    T_sw = np.ones((7, 1)) * setup['Dat']['stat']['Tj']
    T_ca = setup['Dat']['stat']['Tj']
    T_old = setup['Dat']['stat']['Tj']
    iter1 = 0

    # ==============================================================================
    # Variables
    # ==============================================================================
    [_, timeElec, timeLoss, timeTher, _, _, _, _, _] = top.initOut()

    ###################################################################################################################
    # Pre-Processing
    ###################################################################################################################
    # ==============================================================================
    # Generate Reference Waveform
    # ==============================================================================
    t_ref = np.linspace(0, K / fel, K * N + 1)
    [v_ref, e_ref, _] = top.calcRef(E, phiE, phiV, setup)

    # ==============================================================================
    # Thermal ROM
    # ==============================================================================
    # ------------------------------------------
    # Load
    # ------------------------------------------
    [Rth_JA, Cth_JA, Rth_DA, Cth_DA, Rth_CA, Cth_CA, Rth_JA_cap, Cth_JA_cap] = initRC(para, setup)

    # ------------------------------------------
    # Variables
    # ------------------------------------------
    Tinit_T = np.zeros((len(Rth_JA), len(top.id2)))
    Tinit_D = np.zeros((len(Rth_DA), len(top.id2)))
    Tinit_K = np.zeros((len(Rth_CA), len(top.id2)))
    Tinit_C = np.zeros(np.size(Rth_JA_cap))

    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    # ==============================================================================
    # Switching Function
    # ==============================================================================
    [xs, xsh, s, c, x, xN0] = top.calcPWM(v_ref, t_ref, Mi, setup)

    # ==============================================================================
    # Time Domain
    # ==============================================================================
    [timeAc, timeDc, _] = top.calcTime(s, e_ref, t_ref, Mi, mdl, start, ende, [], 1, setup)

    # ==============================================================================
    # Msg
    # ==============================================================================
    print("START: Iterating")
    print("------------------------------------------")

    # ==============================================================================
    # Iterating
    # ==============================================================================
    while err > setup['Exp']['tol']:
        # ------------------------------------------
        # Electrical
        # ------------------------------------------
        # Switches
        for i in range(0, len(top.id2)):
            timeElec['sw'][top.id2[i]] = calcElecSwi(Vdc, top.id9[i] * timeAc[top.id4[i]], (s[top.id3[i]][start:ende] == (-1) ** i),
                                                     T_sw[i], top.id5[i], para, setup)

        # Capacitor
        timeDc['v_dc'] = calcElecCap(t_ref, timeDc['i_c'], T_ca, para, setup)
        timeElec['cap']['C1']['i_c'] = timeDc['i_c']
        timeElec['cap']['C1']['v_c'] = timeDc['v_dc']

        # ------------------------------------------
        # Losses
        # ------------------------------------------
        # Switches
        for i in range(0, len(top.id2)):
            timeLoss['sw'][top.id2[i]] = calcLossSwi(s[top.id3[i]][start:ende] * (-1) ** i,
                                                     timeElec['sw'][top.id2[i]]['i_T'],
                                                     timeElec['sw'][top.id2[i]]['i_D'],
                                                     timeElec['sw'][top.id2[i]]['v_T'],
                                                     timeElec['sw'][top.id2[i]]['v_D'],
                                                     T_sw[i], para, setup)

        # Capacitor
        timeLoss['cap']['C1'] = calcLossCap(t_ref, timeDc['i_c'], T_ca, para, setup)

        # ------------------------------------------
        # Init Thermal
        # ------------------------------------------
        if iter1 == 0:
            # Switches
            for i in range(0, len(top.id2)):
                Tinit_T[:, i] = np.mean(timeLoss['sw'][top.id2[i]]['p_T']) * Rth_JA
                Tinit_D[:, i] = np.mean(timeLoss['sw'][top.id2[i]]['p_D']) * Rth_DA
                if setup['Par']['Ther']['Coupling'] != 0:
                    Tinit_K[:, i] = np.mean(timeLoss['sw'][top.id2[i]]['p_L']) * Rth_CA

            # Capacitor
            Tinit_C = np.mean(timeLoss['cap']['C1']['p_L']) * Rth_JA_cap

        # ------------------------------------------
        # Thermal
        # ------------------------------------------
        # Baseplate
        if setup['Par']['Ther']['Coupling'] == 2:
            Pv_Base = np.mean((timeLoss['sw']['S1']['p_L'] + timeLoss['sw']['S2']['p_L'] + timeLoss['sw']['S3']['p_L'] +
                               timeLoss['sw']['S4']['p_L'] + timeLoss['sw']['S5']['p_L'] + timeLoss['sw']['S6']['p_L'])) / 6
        else:
            Pv_Base = 0

        # Switches
        for i in range(0, len(top.id2)):
            [timeTher['sw'][top.id6[i]], Tinit_T[:, i]] = calcTherRC(Tinit_T[:, i], Tc, timeLoss['sw'][top.id2[i]]['p_T'],
                                                                     t_ref[start:ende], Rth_JA, Cth_JA)
            [timeTher['sw'][top.id7[i]], Tinit_D[:, i]] = calcTherRC(Tinit_D[:, i], Tc, timeLoss['sw'][top.id2[i]]['p_D'],
                                                                     t_ref[start:ende], Rth_DA, Cth_DA)
            if setup['Par']['Ther']['Coupling'] == 1:
                [timeTher['sw'][top.id8[i]], Tinit_K[:, i]] = calcTherRC(Tinit_K[:, i], Tc, timeLoss['sw'][top.id2[i]]['p_L'],
                                                                         t_ref[start:ende], Rth_CA, Cth_CA)
                timeTher['sw'][top.id6[i]] = timeTher['sw'][top.id6[i]][:] + timeTher['sw'][top.id8[i]][:] - Tc
                timeTher['sw'][top.id7[i]] = timeTher['sw'][top.id7[i]][:] + timeTher['sw'][top.id8[i]][:] - Tc
            elif setup['Par']['Ther']['Coupling'] == 2:
                dT_Base = Rth_CA * Pv_Base
                timeTher['sw'][top.id6[i]] = timeTher['sw'][top.id6[i]][:] + dT_Base
                timeTher['sw'][top.id7[i]] = timeTher['sw'][top.id7[i]][:] + dT_Base
                timeTher['sw'][top.id8[i]] = Tc * np.ones(np.size(s['A'][start:ende])) + dT_Base
            else:
                timeTher['sw'][top.id8[i]] = Tc * np.ones(np.size(s['A'][start:ende]))

        # Capacitor
        [timeTher['cap']['C1'], Tinit_C] = calcTherRC(Tinit_C, Tc, timeLoss['cap']['C1']['p_L'], t_ref[start:ende],
                                                      Rth_JA_cap, Cth_JA_cap)

        # ------------------------------------------
        # Error
        # ------------------------------------------
        err = np.mean(np.abs(timeTher['sw']['T1'] - T_old)) / np.mean(T_old)

        # ------------------------------------------
        # Update Para
        # ------------------------------------------
        if setup['Exp']['therFeed'] == 1:
            # Switches
            for i in range(0, len(top.id2)):
                T_sw[i] = np.mean(timeTher['sw'][top.id6[i]])

            # Capacitor
            T_ca = np.mean(timeTher['cap']['C1'])

        # Iteration
        iter1 = iter1 + 1

        # Previous Temperature
        T_old = timeTher['sw']['T1']

        # ------------------------------------------
        # Msg
        # ------------------------------------------
        if iter1 < int(setup['Exp']['int']):
            print(
                "ITER: %d) Stationary temperature T_swi=%.2f C (T_cap=%.2f C) and P_swi=%.2f W (Pv_cap=%.2f W) with error: %.2f %%" % (
                    iter1, T_sw[0], T_ca, np.mean(timeLoss['sw']['S1']['p_L']), np.mean(timeLoss['cap']['C1']['p_L']),
                    err * 100))
        else:
            print("ITER: %d) Maximum iteration reached" % iter1)
            break

    # ==============================================================================
    # Msg
    # ==============================================================================
    print("------------------------------------------")
    print("INFO: Converged after %d iterations with T_swi=%.2f C (T_cap=%.2f C) and error: %.2f %%" % (
        iter1, T_sw[0], T_ca, err * 100))

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
    out = calcAvg(top, out, setup)

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
    # Baseplate
    # ------------------------------------------
    if setup['Par']['Ther']['Coupling'] == 2:
        Pv_Base = np.mean((timeLoss['sw']['S1']['p_L'] + timeLoss['sw']['S2']['p_L'] + timeLoss['sw']['S3']['p_L'] +
                           timeLoss['sw']['S4']['p_L'] + timeLoss['sw']['S5']['p_L'] + timeLoss['sw']['S6']['p_L'])) / 6
    else:
        Pv_Base = 0

    # ------------------------------------------
    # Switches
    # ------------------------------------------
    for i in range(0, len(top.id2)):
        [timeTher['sw'][top.id6[i]], _] = calcTherRC(Tinit_T[:, i], Tc, timeLoss['sw'][top.id2[i]]['p_T'],
                                                     t_ref[start:ende], Rth_JA, Cth_JA)
        [timeTher['sw'][top.id7[i]], _] = calcTherRC(Tinit_D[:, i], Tc, timeLoss['sw'][top.id2[i]]['p_D'],
                                                     t_ref[start:ende], Rth_DA, Cth_DA)
        if setup['Par']['Ther']['Coupling'] == 1:
            [timeTher['sw'][top.id8[i]], _] = calcTherRC(Tinit_K[:, i], Tc, timeLoss['sw'][top.id2[i]]['p_L'],
                                                         t_ref[start:ende], Rth_CA, Cth_CA)
            timeTher['sw'][top.id6[i]] = timeTher['sw'][top.id6[i]][:] + timeTher['sw'][top.id8[i]][:] - Tc
            timeTher['sw'][top.id7[i]] = timeTher['sw'][top.id7[i]][:] + timeTher['sw'][top.id8[i]][:] - Tc
        elif setup['Par']['Ther']['Coupling'] == 2:
            dT_Base = Rth_CA * Pv_Base
            timeTher['sw'][top.id6[i]] = timeTher['sw'][top.id6[i]][:] + dT_Base
            timeTher['sw'][top.id7[i]] = timeTher['sw'][top.id7[i]][:] + dT_Base
            timeTher['sw'][top.id8[i]] = Tc * np.ones(np.size(s['A'][start:ende])) + dT_Base
        else:
            timeTher['sw'][top.id8[i]] = Tc * np.ones(np.size(s['A'][start:ende]))

    # ------------------------------------------
    # Capacitor
    # ------------------------------------------
    [timeTher['cap']['C1'], _] = calcTherRC(Tinit_C, Tc, timeLoss['cap']['C1']['p_L'], t_ref[start:ende], Rth_JA_cap,
                                            Cth_JA_cap)

    # ==============================================================================
    # Frequency domain
    # ==============================================================================
    [freqSw, freqAc, freqDc] = calcFreq(s['A'][start:ende], xs['A'][start:ende], timeAc['i_a'], timeAc['v_a'],
                                        timeAc['v_a0'], timeDc['i_dc'], timeDc['v_dc'])

    # ==============================================================================
    # Output
    # ==============================================================================
    [time, freq, _] = top.out(timeElec, timeLoss, timeTher, timeAc, timeDc, freqSw, freqAc, freqDc, [], [], t_ref, v_ref,
                              e_ref, s, c, xs, xsh, x, xN0, [], start, ende, 1)

    ###################################################################################################################
    # MSG Out
    ###################################################################################################################
    print("------------------------------------------")
    print("END: Steady-State solution class", top.name)
    print("------------------------------------------")

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return [time, freq]
