#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         calcSteady_DCDC
# Date:         08.02.2026
# Author:       Dr. Pascal A. Schirmer
# Version:      V.1.1
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
from src.elec.calcElecTra import calcElecTra
from src.elec.calcLossTra import calcLossTra
from src.general.calcAvg import calcAvg

# ==============================================================================
# External
# ==============================================================================
import numpy as np
import math


#######################################################################################################################
# Function
#######################################################################################################################
def calcSteady_DCDC(top, mdl, para, setup):
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
    T_sw = np.ones((len(top.id2), 1)) * setup['Dat']['stat']['Tj']
    T_ca = np.ones((2, 1)) * setup['Dat']['stat']['Tj']
    T_tra = np.ones((3, 1)) * setup['Dat']['stat']['Tj']
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
    [v_ref, e_ref, _] = top.calcRef(E, phiE, phiV, [], setup)

    # ==============================================================================
    # Thermal ROM
    # ==============================================================================
    # ------------------------------------------
    # Load
    # ------------------------------------------
    # Primary Switches
    para_pri = {'Swi': para['SwiPri'], 'Cap': para['CapPri'], 'Tra': para['Tra']}
    [Rth_JA, Cth_JA, Rth_DA, Cth_DA, Rth_CA, Cth_CA, Rth_JA_cap, Cth_JA_cap,
     Rth_PA_tra, Cth_PA_tra, Rth_SA_tra, Cth_SA_tra, Rth_CA_tra, Cth_CA_tra] = initRC(para_pri, setup)

    # Secondary Switches
    para_sec = {'Swi': para['SwiSec'], 'Cap': para['Cap'], 'Tra': para['Tra']}
    [Rth_JA_s, Cth_JA_s, Rth_DA_s, Cth_DA_s, Rth_CA_s, Cth_CA_s, Rth_JA_cap_s, Cth_JA_cap_s,
     _, _, _, _, _, _] = initRC(para_sec, setup)

    # Sanity Check
    if len(Rth_JA_s) != len(Rth_JA) or len(Rth_DA_s) != len(Rth_DA) or len(Rth_CA_s) != len(Rth_CA):
        print("WARN: DAB primary/secondary RC networks have different lengths; using primary RC lengths for both.")
        Rth_JA_s, Cth_JA_s, Rth_DA_s, Cth_DA_s, Rth_CA_s, Cth_CA_s = Rth_JA, Cth_JA, Rth_DA, Cth_DA, Rth_CA, Cth_CA

    # Sanity Check
    if len(Rth_JA_cap) != len(Rth_JA_cap_s):
        print("WARN: DAB primary/secondary Cap networks have different lengths; using primary RC lengths for both.")
        Rth_JA_cap_s, Cth_JA_s = Cth_JA_cap_s, Cth_JA_cap

    # ------------------------------------------
    # Variables
    # ------------------------------------------
    Tinit_T = np.zeros((len(Rth_JA), len(top.id2)))
    Tinit_D = np.zeros((len(Rth_DA), len(top.id2)))
    Tinit_K = np.zeros((len(Rth_CA), len(top.id2)))
    Tinit_C1 = np.zeros(np.size(Rth_JA_cap))
    Tinit_C2 = np.zeros(np.size(Rth_JA_cap_s))
    Tinit_PC = np.zeros(np.size(Rth_PA_tra))
    Tinit_SC = np.zeros(np.size(Rth_SA_tra))
    Tinit_CC = np.zeros(np.size(Rth_CA_tra))

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
    [timeAc, timeDc, _] = top.calcTime(s, e_ref, t_ref, Mi, mdl, start, ende, [], para, setup)

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
            para_swi = para
            if 'SwiPri' in para and 'SwiSec' in para:
                para_swi = {'Swi': para['SwiPri'], 'Cap': para['Cap']} if i < 4 else {'Swi': para['SwiSec'], 'Cap': para['Cap']}

            timeElec['sw'][top.id2[i]] = calcElecSwi(Vdc, top.id9[i] * timeAc[top.id4[i]], (s[top.id3[i]][start:ende] == (-1) ** i),
                                                     T_sw[i], top.id5[i], para_swi, setup)

        # Output Capacitor
        timeElec['cap']['C2']['i_c'] = timeDc['i_c_sec']
        timeElec['cap']['C2']['v_c'] = timeDc['v_dc_sec_cap']

        # Input Capacitor
        timeElec['cap']['C1']['i_c'] = timeDc['i_c_pri']
        timeElec['cap']['C1']['v_c'] = timeDc['v_dc_pri_cap']

        # Transformer
        timeElec['tra'] = calcElecTra(t_ref[start:ende], timeAc['i_ac_pri'], timeAc['i_ac_sec'], timeAc['v_ac_pri'], timeAc['v_ac_sec'], T_tra, para, setup)

        # ------------------------------------------
        # Losses
        # ------------------------------------------
        # Switches
        for i in range(0, len(top.id2)):
            para_swi = para
            if 'SwiPri' in para and 'SwiSec' in para:
                para_swi = {'Swi': para['SwiPri'], 'Cap': para['CapPri']} if i < 4 else {'Swi': para['SwiSec'], 'Cap': para['CapSec']}

            timeLoss['sw'][top.id2[i]] = calcLossSwi(s[top.id3[i]][start:ende] * (-1) ** i,
                                                     timeElec['sw'][top.id2[i]]['i_T'],
                                                     timeElec['sw'][top.id2[i]]['i_D'],
                                                     timeElec['sw'][top.id2[i]]['v_T'],
                                                     timeElec['sw'][top.id2[i]]['v_D'],
                                                     T_sw[i], para_swi, setup)

        # Capacitor
        timeLoss['cap']['C1'] = calcLossCap(t_ref, timeElec['cap']['C1']['i_c'], T_ca[0], para, setup)
        timeLoss['cap']['C2'] = calcLossCap(t_ref, timeElec['cap']['C2']['i_c'], T_ca[1], para, setup)

        # Transformer
        timeLoss['tra'] = calcLossTra(t_ref[start:ende], timeElec['tra']['i_p'], timeElec['tra']['i_s'],
                                      timeElec['tra']['v_p'], timeElec['tra']['v_s'], T_tra, para, setup)

        # ------------------------------------------
        # Init Thermal
        # ------------------------------------------
        if iter1 == 0:
            # Switches
            for i in range(0, len(top.id2)):
                if i >= 4:
                    Tinit_T[:, i] = np.mean(timeLoss['sw'][top.id2[i]]['p_T']) * Rth_JA_s
                    Tinit_D[:, i] = np.mean(timeLoss['sw'][top.id2[i]]['p_D']) * Rth_DA_s
                    if setup['Par']['Ther']['Coupling'] != 0:
                        Tinit_K[:, i] = np.mean(timeLoss['sw'][top.id2[i]]['p_L']) * Rth_CA_s
                else:
                    Tinit_T[:, i] = np.mean(timeLoss['sw'][top.id2[i]]['p_T']) * Rth_JA
                    Tinit_D[:, i] = np.mean(timeLoss['sw'][top.id2[i]]['p_D']) * Rth_DA
                    if setup['Par']['Ther']['Coupling'] != 0:
                        Tinit_K[:, i] = np.mean(timeLoss['sw'][top.id2[i]]['p_L']) * Rth_CA

            # Capacitor
            Tinit_C1 = np.mean(timeLoss['cap']['C1']['p_L']) * Rth_JA_cap
            Tinit_C2 = np.mean(timeLoss['cap']['C2']['p_L']) * Rth_JA_cap_s

            # Transformer
            Tinit_PC = np.mean(timeLoss['tra']['p_PC']) * Rth_PA_tra
            Tinit_SC = np.mean(timeLoss['tra']['p_SC']) * Rth_SA_tra
            Tinit_CC = np.mean(timeLoss['tra']['p_CC']) * Rth_CA_tra

        # ------------------------------------------
        # Thermal
        # ------------------------------------------
        # Baseplate
        if setup['Par']['Ther']['Coupling'] == 2:
            Pv_Base_all = np.zeros((2, 1))
            for i in range(0, len(top.id2)):
                if i < 4:
                    Pv_Base_all[0] = Pv_Base_all[0] + np.mean(timeLoss['sw'][top.id2[i]]['p_L'])
                else:
                    Pv_Base_all[1] = Pv_Base_all[1] + np.mean(timeLoss['sw'][top.id2[i]]['p_L'])
            Pv_Base_all = Pv_Base_all / (len(top.id2) / 2)
        else:
            Pv_Base_all = np.zeros((2, 1))

        # Switches
        for i in range(0, len(top.id2)):
            if i >= 4:
                rth_ja, cth_ja, rth_da, cth_da, rth_ca, cth_ca = Rth_JA_s, Cth_JA_s, Rth_DA_s, Cth_DA_s, Rth_CA_s, Cth_CA_s
                Pv_Base = Pv_Base_all[1]
            else:
                rth_ja, cth_ja, rth_da, cth_da, rth_ca, cth_ca = Rth_JA, Cth_JA, Rth_DA, Cth_DA, Rth_CA, Cth_CA
                Pv_Base = Pv_Base_all[0]

            [timeTher['sw'][top.id6[i]], Tinit_T[:, i]] = calcTherRC(Tinit_T[:, i], Tc, timeLoss['sw'][top.id2[i]]['p_T'],
                                                                     t_ref[start:ende], rth_ja, cth_ja)
            [timeTher['sw'][top.id7[i]], Tinit_D[:, i]] = calcTherRC(Tinit_D[:, i], Tc, timeLoss['sw'][top.id2[i]]['p_D'],
                                                                     t_ref[start:ende], rth_da, cth_da)
            if setup['Par']['Ther']['Coupling'] == 1:
                [timeTher['sw'][top.id8[i]], Tinit_K[:, i]] = calcTherRC(Tinit_K[:, i], Tc, timeLoss['sw'][top.id2[i]]['p_L'],
                                                                         t_ref[start:ende], rth_ca, cth_ca)
                timeTher['sw'][top.id6[i]] = timeTher['sw'][top.id6[i]][:] + timeTher['sw'][top.id8[i]][:] - Tc
                timeTher['sw'][top.id7[i]] = timeTher['sw'][top.id7[i]][:] + timeTher['sw'][top.id8[i]][:] - Tc
            elif setup['Par']['Ther']['Coupling'] == 2:
                dT_Base = rth_ca * Pv_Base
                timeTher['sw'][top.id6[i]] = timeTher['sw'][top.id6[i]][:] + dT_Base
                timeTher['sw'][top.id7[i]] = timeTher['sw'][top.id7[i]][:] + dT_Base
                timeTher['sw'][top.id8[i]] = Tc * np.ones(np.size(s['A'][start:ende])) + dT_Base
            else:
                timeTher['sw'][top.id8[i]] = Tc * np.ones(np.size(s['A'][start:ende]))

        # Capacitor
        [timeTher['cap']['C1'], Tinit_C1] = calcTherRC(Tinit_C1, Tc, timeLoss['cap']['C1']['p_L'], t_ref[start:ende],
                                                       Rth_JA_cap, Cth_JA_cap)
        [timeTher['cap']['C2'], Tinit_C2] = calcTherRC(Tinit_C2, Tc, timeLoss['cap']['C2']['p_L'], t_ref[start:ende],
                                                       Rth_JA_cap_s, Cth_JA_cap_s)

        # Transformer
        [timeTher['tra']['PC'], Tinit_PC] = calcTherRC(Tinit_PC, Tc, timeLoss['tra']['p_PC'], t_ref[start:ende],
                                                       Rth_PA_tra, Cth_PA_tra)
        [timeTher['tra']['SC'], Tinit_SC] = calcTherRC(Tinit_SC, Tc, timeLoss['tra']['p_SC'], t_ref[start:ende],
                                                       Rth_SA_tra, Cth_SA_tra)
        [timeTher['tra']['CC'], Tinit_CC] = calcTherRC(Tinit_CC, Tc, timeLoss['tra']['p_CC'], t_ref[start:ende],
                                                       Rth_CA_tra, Cth_CA_tra)

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
            T_ca[0] = np.mean(timeTher['cap']['C1'])
            T_ca[1] = np.mean(timeTher['cap']['C2'])

            # Transformer
            T_tra[0] = np.mean(timeTher['tra']['PC'])
            T_tra[1] = np.mean(timeTher['tra']['SC'])
            T_tra[2] = np.mean(timeTher['tra']['CC'])

        # Iteration
        iter1 = iter1 + 1

        # Previous Temperature
        T_old = timeTher['sw']['T1']

        # ------------------------------------------
        # Msg
        # ------------------------------------------
        if iter1 < int(setup['Exp']['int']):
            print(
                "ITER: %d) Stationary temperature T_pri=%.2f C (T_sec=%.2f C) and P_pri=%.2f W (Pv_sec=%.2f W) with error: %.2f %%" % (
                    iter1, T_sw[0], T_sw[4], np.mean(timeLoss['sw']['S1']['p_L']), np.mean(timeLoss['sw']['S5']['p_L']),
                    err * 100))
        else:
            print("ITER: %d) Maximum iteration reached" % iter1)
            break

    # ==============================================================================
    # Msg
    # ==============================================================================
    print("------------------------------------------")
    print("INFO: Converged after %d iterations with T_pri=%.2f C (T_sec=%.2f C) and error: %.2f %%" % (
        iter1, T_sw[0], T_sw[4], err * 100))

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
        Pv_Base_all = np.zeros((2, 1))
        for i in range(0, len(top.id2)):
            if i < 4:
                Pv_Base_all[0] = Pv_Base_all[0] + np.mean(timeLoss['sw'][top.id2[i]]['p_L'])
            else:
                Pv_Base_all[1] = Pv_Base_all[1] + np.mean(timeLoss['sw'][top.id2[i]]['p_L'])
        Pv_Base_all = Pv_Base_all / (len(top.id2) / 2)
    else:
        Pv_Base_all = np.zeros((2, 1))

    # ------------------------------------------
    # Switches
    # ------------------------------------------
    for i in range(0, len(top.id2)):
        # Init RC Pairs
        if i >= 4:
            rth_ja, cth_ja, rth_da, cth_da, rth_ca, cth_ca = Rth_JA_s, Cth_JA_s, Rth_DA_s, Cth_DA_s, Rth_CA_s, Cth_CA_s
            Pv_Base = Pv_Base_all[1]
        else:
            rth_ja, cth_ja, rth_da, cth_da, rth_ca, cth_ca = Rth_JA, Cth_JA, Rth_DA, Cth_DA, Rth_CA, Cth_CA
            Pv_Base = Pv_Base_all[0]

        # Calc Transistor and Diode Temperatures
        [timeTher['sw'][top.id6[i]], _] = calcTherRC(Tinit_T[:, i], Tc, timeLoss['sw'][top.id2[i]]['p_T'],
                                                     t_ref[start:ende], rth_ja, cth_ja)
        [timeTher['sw'][top.id7[i]], _] = calcTherRC(Tinit_D[:, i], Tc, timeLoss['sw'][top.id2[i]]['p_D'],
                                                     t_ref[start:ende], rth_da, cth_da)

        # Calculate coupled temperature
        if setup['Par']['Ther']['Coupling'] == 1:
            [timeTher['sw'][top.id8[i]], _] = calcTherRC(Tinit_K[:, i], Tc, timeLoss['sw'][top.id2[i]]['p_L'],
                                                         t_ref[start:ende], rth_ca, cth_ca)
            timeTher['sw'][top.id6[i]] = timeTher['sw'][top.id6[i]][:] + timeTher['sw'][top.id8[i]][:] - Tc
            timeTher['sw'][top.id7[i]] = timeTher['sw'][top.id7[i]][:] + timeTher['sw'][top.id8[i]][:] - Tc
        elif setup['Par']['Ther']['Coupling'] == 2:
            dT_Base = rth_ca * Pv_Base
            timeTher['sw'][top.id6[i]] = timeTher['sw'][top.id6[i]][:] + dT_Base
            timeTher['sw'][top.id7[i]] = timeTher['sw'][top.id7[i]][:] + dT_Base
            timeTher['sw'][top.id8[i]] = Tc * np.ones(np.size(s['A'][start:ende])) + dT_Base
        else:
            timeTher['sw'][top.id8[i]] = Tc * np.ones(np.size(s['A'][start:ende]))

    # ------------------------------------------
    # Capacitor
    # ------------------------------------------
    [timeTher['cap']['C1'], _] = calcTherRC(Tinit_C1, Tc, timeLoss['cap']['C1']['p_L'], t_ref[start:ende], Rth_JA_cap,
                                            Cth_JA_cap)
    [timeTher['cap']['C2'], _] = calcTherRC(Tinit_C2, Tc, timeLoss['cap']['C2']['p_L'], t_ref[start:ende], Rth_JA_cap_s,
                                            Cth_JA_cap_s)

    # ------------------------------------------
    # Transformer
    # ------------------------------------------
    [timeTher['tra']['PC'], _] = calcTherRC(Tinit_PC, Tc, timeLoss['tra']['p_PC'], t_ref[start:ende], Rth_PA_tra, Cth_PA_tra)
    [timeTher['tra']['SC'], _] = calcTherRC(Tinit_SC, Tc, timeLoss['tra']['p_SC'], t_ref[start:ende], Rth_SA_tra, Cth_SA_tra)
    [timeTher['tra']['CC'], _] = calcTherRC(Tinit_CC, Tc, timeLoss['tra']['p_CC'], t_ref[start:ende], Rth_CA_tra, Cth_CA_tra)

    # ==============================================================================
    # Frequency domain
    # ==============================================================================
    [freqSw, freqAc, freqDc] = calcFreq(s['A'][start:ende], xs['A'][start:ende], timeAc['i_ac_pri'], timeAc['v_ac_pri'],
                                        timeAc['v_ac_sec'], timeDc['i_dc_sec'], timeDc['v_dc_sec'])

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
