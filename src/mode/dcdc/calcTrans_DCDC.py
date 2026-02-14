#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         calcTrans_DCDC
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
This function calculates the results for a transient solution for any given topology class.
Inputs:     1) top:     topology class
            2) mdl:     all models and transfer functions of the architecture
            3) para:    all parameters used in the simulation
            4) setup:   includes all simulation variables
Outputs:    1) time:    results in the time domain
            2) freq:    results in the frequency domain
"""

#######################################################################################################################
# Import libs
#######################################################################################################################
# ==============================================================================
# Internal
# ==============================================================================
from src.general.calcSpec import calcFreq
from src.elec.calcElecSwi import calcElecSwi
from src.elec.calcLossSwi import calcLossSwi
from src.therm.calcTherRC import calcTherRC
from src.elec.calcLossCap import calcLossCap
from src.general.calcAvg import calcAvg
from src.therm.initRC import initRC
from src.elec.calcElecCap import calcElecCap
from src.general.append import app_fel, app_fs
from src.elec.calcElecTra import calcElecTra
from src.elec.calcLossTra import calcLossTra

# ==============================================================================
# External
# ==============================================================================
import numpy as np
import math
import pandas as pd
from tqdm import tqdm


#######################################################################################################################
# Function
#######################################################################################################################
def calcTrans_DCDC(top, mdl, para, setup):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("------------------------------------------")
    print("START: Transient solution class", top.name)
    print("------------------------------------------")

    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Parameters
    # ==============================================================================
    fel = setup['Top']['fel']
    fsim = setup['Exp']['fsim']
    fs = setup['Par']['PWM']['fs']
    Tel = 1 / fel
    Nsim = int(np.ceil(fsim / fel))
    Npwm = int(np.ceil(fs / fel))
    K = int(setup['Dat']['stat']['cyc'])
    Nel = int(np.ceil(setup['Dat']['trans']['tmax'] * fel))
    Mi = setup['Dat']['stat']['Mi']

    # ==============================================================================
    # Variables
    # ==============================================================================
    # ------------------------------------------
    # Electrical
    # ------------------------------------------
    E = setup['Top']['E']
    Vdc = setup['Dat']['stat']['Vdc']
    phiE = math.radians(setup['Top']['phiE'])
    phiV = math.radians(setup['Dat']['stat']['phi'])
    Ntr = top.Ntr

    # ------------------------------------------
    # Thermal
    # ------------------------------------------
    Ta = setup['Dat']['stat']['Tc']
    Tj = setup['Dat']['stat']['Tj'] * np.ones((len(top.id2), 1))
    Tcap = setup['Dat']['stat']['Tj']
    T_tra = np.ones((3, 1)) * setup['Dat']['stat']['Tj']

    # ==============================================================================
    # Update Frequency
    # ==============================================================================
    if setup['Exp']['therFeed'] == 0:
        iterNel = 1
        iterNpwm = 1
    elif setup['Exp']['freqPar'] == 'fel':
        iterNel = Nel
        iterNpwm = 1
    else:
        iterNel = Nel
        iterNpwm = Npwm

    # ==============================================================================
    # Outputs
    # ==============================================================================
    out = top.initData()

    ###################################################################################################################
    # Pre-Processing
    ###################################################################################################################
    # ==============================================================================
    # Generate Reference Waveform
    # ==============================================================================
    # ------------------------------------------
    # Time
    # ------------------------------------------
    t_ref = np.linspace(0, K / fel, K * Nsim + 1)

    # ------------------------------------------
    # Reference
    # ------------------------------------------
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
    # Tinit_D = np.zeros((len(Rth_DA), len(top.id2)))
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
    # Ideal Primary and Secondary Voltages
    # ==============================================================================
    v_ac_pri = 0.5 * Vdc * (s['A'] - s['B'])
    v_ac_sec = 0.5 * (Vdc / Ntr) * (s['C'] - s['D']) if Ntr != 0 else 0.5 * Vdc * (s['C'] - s['D'])

    # ==============================================================================
    # Time Domain
    # ==============================================================================
    [timeAc, timeDc, initL] = top.calcTime(s, v_ac_pri, v_ac_sec, e_ref, t_ref, Mi, mdl, Nsim * (K - 1), (K * Nsim + 1), [], para, setup)

    # ==============================================================================
    # Electrical cycle
    # ==============================================================================
    for _ in tqdm(range(iterNel), desc='Elec-Period', position=0):
        # ------------------------------------------
        # Init
        # ------------------------------------------
        dataFel = top.initData()

        # ------------------------------------------
        # Electrical
        # ------------------------------------------
        [timeAc, timeDc, initL] = top.calcTime(s, v_ac_pri, v_ac_sec, e_ref, t_ref, Mi, mdl, Nsim * (K - 1), (K * Nsim + 1), initL, para, setup)
        # timeDc['v_dc'] = calcElecCap(t_ref, timeDc['i_c'], Tcap, para, setup)

        # ------------------------------------------
        # PWM Period
        # ------------------------------------------
        for ii in range(iterNpwm):
            # Init
            [_, timeElec, timeLoss, timeTher, _, _, _, _, _] = top.initOut()
            start = int(ii * (Nsim / iterNpwm))
            if iterNpwm == 1:
                ende = int(Nsim / iterNpwm * (ii + 1) + 1)
            else:
                ende = int(Nsim / iterNpwm * (ii + 1) + 0)

            # Switch
            for j in range(0, len(top.id2)):
                para_swi = para
                Vdc_temp = Vdc
                if 'SwiPri' in para and 'SwiSec' in para:
                    para_swi = {'Swi': para['SwiPri'], 'Cap': para.get('CapPri', para.get('Cap', []))} if j < 4 else {'Swi': para['SwiSec'], 'Cap': para.get('CapSec', para.get('Cap', []))}
                    Vdc_temp = np.mean(timeDc['v_dc_pri']) if j < 4 else np.mean(timeDc['v_dc_sec'])

                timeElec['sw'][top.id2[j]] = calcElecSwi(Vdc_temp, top.id9[j] * timeAc[top.id4[j]][start:ende], (s[top.id3[j]][start:ende] == (-1) ** j), Tj[j], top.id5[j], para_swi, setup)
                timeLoss['sw'][top.id2[j]] = calcLossSwi(s[top.id3[j]][start:ende] * (-1) ** j, timeElec['sw'][top.id2[j]]['i_T'], timeElec['sw'][top.id2[j]]['i_D'], timeElec['sw'][top.id2[j]]['v_T'], timeElec['sw'][top.id2[j]]['v_D'], Tj[j], para_swi, setup)

                if setup['Par']['Ther']['Heatsink'] == 1 and setup['Par']['Ther']['Coupling'] == 1:
                    if j >= 4:
                        rth_ja, cth_ja, rth_ca, cth_ca = Rth_JA_s, Cth_JA_s, Rth_CA_s, Cth_CA_s
                    else:
                        rth_ja, cth_ja, rth_ca, cth_ca = Rth_JA, Cth_JA, Rth_CA, Cth_CA
                    [timeTher['sw'][top.id6[j]], Tinit_T[:, j]] = calcTherRC(Tinit_T[:, j], Ta, timeLoss['sw'][top.id2[j]]['p_T'], t_ref[start:ende], rth_ja, cth_ja)
                    [timeTher['sw'][top.id8[j]], Tinit_K[:, j]] = calcTherRC(Tinit_K[:, j], Ta, timeLoss['sw'][top.id2[j]]['p_L'], t_ref[start:ende], rth_ca, cth_ca)
                    timeTher['sw'][top.id6[j]] = timeTher['sw'][top.id6[j]][:] + timeTher['sw'][top.id8[j]][:] - Ta
                elif setup['Par']['Ther']['Heatsink'] == 1 and setup['Par']['Ther']['Coupling'] == 2:
                    if j >= 4:
                        rth_ja, cth_ja, rth_ca, cth_ca = Rth_JA_s, Cth_JA_s, Rth_CA_s, Cth_CA_s
                    else:
                        rth_ja, cth_ja, rth_ca, cth_ca = Rth_JA, Cth_JA, Rth_CA, Cth_CA
                    [timeTher['sw'][top.id6[j]], Tinit_T[:, j]] = calcTherRC(Tinit_T[:, j], Ta, timeLoss['sw'][top.id2[j]]['p_T'], t_ref[start:ende], rth_ja, cth_ja)
                    [timeTher['sw'][top.id8[j]], Tinit_K[:, j]] = calcTherRC(Tinit_K[:, j], Ta, np.mean(timeLoss['sw'][top.id2[j]]['p_L'])*np.ones(np.size(timeLoss['sw'][top.id2[j]]['p_L'])), t_ref[start:ende], rth_ca, cth_ca)
                    timeTher['sw'][top.id6[j]] = timeTher['sw'][top.id6[j]][:] + timeTher['sw'][top.id8[j]][:] - Ta
                else:
                    if j >= 4:
                        rth_ja, cth_ja = Rth_JA_s, Cth_JA_s
                    else:
                        rth_ja, cth_ja = Rth_JA, Cth_JA
                    [timeTher['sw'][top.id6[j]], Tinit_T[:, j]] = calcTherRC(Tinit_T[:, j], Ta, timeLoss['sw'][top.id2[j]]['p_T'], t_ref[start:ende], rth_ja, cth_ja)

            # Capacitor Primary capacitor
            timeElec['cap']['C1']['i_c'] = timeDc['i_c_pri'][start:ende]
            timeElec['cap']['C1']['v_c'] = timeDc['v_dc_pri'][start:ende]
            timeLoss['cap']['C1'] = calcLossCap(t_ref[start:ende], timeDc['i_c_pri'][start:ende], Tcap, {'Cap': para.get('CapPri', para.get('Cap', []))}, setup)
            [timeTher['cap']['C1'], Tinit_C1] = calcTherRC(Tinit_C1, Ta, timeLoss['cap']['C1']['p_L'], t_ref[start:ende], Rth_JA_cap, Cth_JA_cap)

            # Capacitor Secondary capacitor
            timeElec['cap']['C2']['i_c'] = timeDc['i_c_sec'][start:ende]
            timeElec['cap']['C2']['v_c'] = timeDc['v_dc_sec'][start:ende]
            timeLoss['cap']['C2'] = calcLossCap(t_ref[start:ende], timeDc['i_c_sec'][start:ende], Tcap, {'Cap': para.get('CapSec', para.get('Cap', []))}, setup)
            [timeTher['cap']['C2'], Tinit_C2] = calcTherRC(Tinit_C2, Ta, timeLoss['cap']['C2']['p_L'], t_ref[start:ende], Rth_JA_cap_s, Cth_JA_cap_s)

            # Transformer
            timeElec['tra']['T1'] = calcElecTra(t_ref[start:ende], timeAc['i_ac_pri'][start:ende], timeAc['i_ac_sec'][start:ende],
                                                timeAc['v_ac_pri'][start:ende], timeAc['v_ac_sec'][start:ende], T_tra, para, setup)
            timeLoss['tra']['T1'] = calcLossTra(t_ref[start:ende], timeElec['tra']['T1']['i_p'], timeElec['tra']['T1']['i_s'],
                                                timeElec['tra']['T1']['v_p'], timeElec['tra']['T1']['v_s'], T_tra, para, setup)
            [timeTher['tra']['PC'], Tinit_PC] = calcTherRC(Tinit_PC, Ta, timeLoss['tra']['T1']['p_PC'], t_ref[start:ende],
                                                           Rth_PA_tra, Cth_PA_tra)
            [timeTher['tra']['SC'], Tinit_SC] = calcTherRC(Tinit_SC, Ta, timeLoss['tra']['T1']['p_SC'], t_ref[start:ende],
                                                           Rth_SA_tra, Cth_SA_tra)
            [timeTher['tra']['CC'], Tinit_CC] = calcTherRC(Tinit_CC, Ta, timeLoss['tra']['T1']['p_CC'], t_ref[start:ende],
                                                           Rth_CA_tra, Cth_CA_tra)

            # Appending
            dataFel = app_fs(dataFel, timeElec, timeLoss, setup)

            # Parameter Update
            if setup['Exp']['therFeed'] == 1:
                for j in range(0, len(top.id2)):
                    Tj[j] = timeTher['sw'][top.id6[j]][-1]
                Tcap = timeTher['cap']['C1'][-1]

        # ------------------------------------------
        # Appending
        # ------------------------------------------
        out = app_fel(out, dataFel['elec'], dataFel['loss'], Nel, setup)

        # ------------------------------------------
        # Update AC Voltage
        # ------------------------------------------
        v_ac_pri[Nsim * (K - 1):(K * Nsim + 1)] = (dataFel['elec']['sw']['S2']['v_T'] - dataFel['elec']['sw']['S1']['v_T'])
        v_ac_sec[Nsim * (K - 1):(K * Nsim + 1)] = (dataFel['elec']['sw']['S6']['v_T'] - dataFel['elec']['sw']['S5']['v_T'])
        timeAc['v_ac_pri'] = v_ac_pri
        timeAc['v_ac_sec'] = v_ac_sec

    # ==============================================================================
    # Averaging
    # ==============================================================================
    out = calcAvg(top, out, setup)

    # ==============================================================================
    # Thermal
    # ==============================================================================
    # ------------------------------------------
    # Init
    # ------------------------------------------
    t = np.linspace(0, Tel * Nel, len(out['loss']['sw'][top.id2[0]]['p_T']))

    # ------------------------------------------
    # Calc
    # ------------------------------------------
    # Switches
    for i in range(0, len(top.id2)):
        if i >= 4:
            rth_ja, cth_ja, rth_da, cth_da, rth_ca, cth_ca = Rth_JA_s, Cth_JA_s, Rth_DA_s, Cth_DA_s, Rth_CA_s, Cth_CA_s
        else:
            rth_ja, cth_ja, rth_da, cth_da, rth_ca, cth_ca = Rth_JA, Cth_JA, Rth_DA, Cth_DA, Rth_CA, Cth_CA

        [out['ther']['sw'][top.id6[i]], _] = calcTherRC(0, Ta, out['loss']['sw'][top.id2[i]]['p_T'].values, t, rth_ja, cth_ja)
        [out['ther']['sw'][top.id7[i]], _] = calcTherRC(0, Ta, out['loss']['sw'][top.id2[i]]['p_D'].values, t, rth_da, cth_da)

    # Capacitor
    [out['ther']['cap']['C1'], _] = calcTherRC(0, Ta, out['loss']['cap']['C1']['p_L'].values, t, Rth_JA_cap, Cth_JA_cap)
    [out['ther']['cap']['C2'], _] = calcTherRC(0, Ta, out['loss']['cap']['C2']['p_L'].values, t, Rth_JA_cap_s, Cth_JA_cap_s)

    # Transformer
    [out['ther']['tra']['PC'], _] = calcTherRC(0, Ta, out['loss']['tra']['T1']['p_PC'].values, t, Rth_PA_tra, Cth_PA_tra)
    [out['ther']['tra']['SC'], _] = calcTherRC(0, Ta, out['loss']['tra']['T1']['p_SC'].values, t, Rth_SA_tra, Cth_SA_tra)
    [out['ther']['tra']['CC'], _] = calcTherRC(0, Ta, out['loss']['tra']['T1']['p_CC'].values, t, Rth_CA_tra, Cth_CA_tra)

    # Coupling
    for i in range(0, len(top.id2)):
        if setup['Par']['Ther']['Heatsink'] == 1 and setup['Par']['Ther']['Coupling'] == 1:
            [out['ther']['sw'][top.id8[i]], _] = calcTherRC(0, Ta, out['loss']['sw'][top.id2[i]]['p_L'].values, t, rth_ca, cth_ca)
            out['ther']['sw'][top.id6[i]] = out['ther']['sw'][top.id6[i]][:] + out['ther']['sw'][top.id8[i]][:] - Ta
            out['ther']['sw'][top.id7[i]] = out['ther']['sw'][top.id7[i]][:] + out['ther']['sw'][top.id8[i]][:] - Ta
        elif setup['Par']['Ther']['Heatsink'] == 1 and setup['Par']['Ther']['Coupling'] == 2:
            Pv_Base = 0
            for sid in top.id2:
                Pv_Base = Pv_Base + out['loss']['sw'][sid]['p_L']
            Pv_Base = np.mean(Pv_Base) / len(top.id2)
            [out['ther']['sw'][top.id8[i]], _] = calcTherRC(0, Ta, Pv_Base*np.ones(np.size(out['loss']['sw'][top.id2[0]]['p_L'])), t, rth_ca, cth_ca)
            out['ther']['sw'][top.id6[i]] = out['ther']['sw'][top.id6[i]][:] + out['ther']['sw'][top.id8[i]][:] - Ta
            out['ther']['sw'][top.id7[i]] = out['ther']['sw'][top.id7[i]][:] + out['ther']['sw'][top.id8[i]][:] - Ta
        else:
            out['ther']['sw'][top.id8[i]] = Ta

    # ------------------------------------------
    # Appending
    # ------------------------------------------
    # Switch
    out['ther']['sw'] = pd.DataFrame(out['ther']['sw'], columns=top.id6 + top.id7 + top.id8)

    # Capacitor
    out['ther']['cap'] = pd.DataFrame(out['ther']['cap'], columns=['C1', 'C2'])

    # Transformer
    out['ther']['tra'] = pd.DataFrame(out['ther']['tra'], columns=['PC', 'SC', 'CC'])

    ###################################################################################################################
    # Post-Processing
    ###################################################################################################################
    # ==============================================================================
    # Frequency domain
    # ==============================================================================
    # ------------------------------------------
    # Calculation
    # ------------------------------------------
    dictSw = {'Sa': s['A'][Nsim:(K * Nsim + 1)], 'Xas': xs['A'][Nsim:(K * Nsim + 1)]}
    if setup['Top']['sourceType'] == 'DAB':
        dictSw['Sb'] = s['C'][Nsim:(K * Nsim + 1)]
    [freqSw, freqAc, freqDc] = calcFreq(dictSw, timeAc, timeDc)

    # ==============================================================================
    # Output
    # ==============================================================================
    [time, freq, _] = top.out(out['elec'], out['loss'], out['ther'], timeAc, timeDc, freqSw, freqAc, freqDc, [], [], t_ref,
                              v_ref, e_ref, s, c, xs, xsh, x, xN0, [], Nsim * (K - 1), (K * Nsim + 1), 1)

    # Re-align lengths because top.out might have different internal logic for slicing
    time['Loss'] = out['loss']
    time['Elec'] = out['elec']
    time['Ther'] = out['ther']

    ###################################################################################################################
    # MSG Out
    ###################################################################################################################
    print("------------------------------------------")
    print("END: Transient class", top.name)
    print("------------------------------------------")

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return [time, freq]
