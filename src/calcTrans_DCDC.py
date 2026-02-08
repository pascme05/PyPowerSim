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
from src.general.calcFreq import calcFreq
from src.elec.calcElecSwi import calcElecSwi
from src.elec.calcLossSwi import calcLossSwi
from src.therm.calcTherRC import calcTherRC
from src.elec.calcLossCap import calcLossCap
from src.general.calcAvg import calcAvg
from src.therm.initRC import initRC
from src.elec.calcElecCap import calcElecCap
from src.general.append import app_fel, app_fs

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

    # ------------------------------------------
    # Thermal
    # ------------------------------------------
    Ta = setup['Dat']['stat']['Tc']
    Tj = setup['Dat']['stat']['Tj'] * np.ones((len(top.id2), 1))
    Tcap = setup['Dat']['stat']['Tj']

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
    # Parameters
    # ------------------------------------------
    [Rth_JA, Cth_JA, Rth_DA, Cth_DA, Rth_CA, Cth_CA, Rth_JA_cap, Cth_JA_cap,
     Rth_PC_tra, Cth_PC_tra, Rth_SC_tra, Cth_SC_tra, Rth_CC_tra, Cth_CC_tra] = initRC(para, setup)

    # Secondary switch thermal parameters (DAB only)
    if setup['Top']['sourceType'] == 'DAB' and 'SwiSec' in para:
        para_sec = {'Swi': para['SwiSec'], 'Cap': para['Cap'], 'Tra': para.get('Tra', [])}
        [Rth_JA_s, Cth_JA_s, Rth_DA_s, Cth_DA_s, Rth_CA_s, Cth_CA_s, _, _, _, _, _, _, _, _] = initRC(para_sec, setup)
    else:
        Rth_JA_s, Cth_JA_s, Rth_DA_s, Cth_DA_s, Rth_CA_s, Cth_CA_s = Rth_JA, Cth_JA, Rth_DA, Cth_DA, Rth_CA, Cth_CA

    if len(Rth_JA_s) != len(Rth_JA) or len(Rth_DA_s) != len(Rth_DA) or len(Rth_CA_s) != len(Rth_CA):
        print("WARN: DAB primary/secondary RC networks have different lengths; using primary RC lengths for both.")
        Rth_JA_s, Cth_JA_s, Rth_DA_s, Cth_DA_s, Rth_CA_s, Cth_CA_s = Rth_JA, Cth_JA, Rth_DA, Cth_DA, Rth_CA, Cth_CA

    # ------------------------------------------
    # Init
    # ------------------------------------------
    Tinit_T = np.zeros((len(Rth_JA), len(top.id2)))
    Tinit_C = np.zeros((len(Rth_CA), len(top.id2)))
    Tinit_Cap = np.zeros(np.size(Rth_JA_cap))

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
    [timeAc, timeDc, _] = top.calcTime(s, e_ref, t_ref, Mi, mdl, Nsim * (K - 1), (K * Nsim + 1), [], 1, setup)

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
        timeDc['v_dc'] = calcElecCap(t_ref, timeDc['i_c'], Tcap, para, setup)

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
                if setup['Top']['sourceType'] == 'DAB' and 'SwiPri' in para and 'SwiSec' in para:
                    para_swi = {'Swi': para['SwiPri'], 'Cap': para['Cap']} if j < 4 else {'Swi': para['SwiSec'], 'Cap': para['Cap']}

                timeElec['sw'][top.id2[j]] = calcElecSwi(Vdc, top.id9[j] * timeAc[top.id4[j]][start:ende], (s[top.id3[j]][start:ende] == (-1) ** j), Tj[j], top.id5[j], para_swi, setup)
                timeLoss['sw'][top.id2[j]] = calcLossSwi(s[top.id3[j]][start:ende] * (-1) ** j, timeElec['sw'][top.id2[j]]['i_T'], timeElec['sw'][top.id2[j]]['i_D'], timeElec['sw'][top.id2[j]]['v_T'], timeElec['sw'][top.id2[j]]['v_D'], Tj[j], para_swi, setup)

                if setup['Par']['Ther']['Heatsink'] == 1 and setup['Par']['Ther']['Coupling'] == 1:
                    if setup['Top']['sourceType'] == 'DAB' and j >= 4:
                        rth_ja, cth_ja, rth_ca, cth_ca = Rth_JA_s, Cth_JA_s, Rth_CA_s, Cth_CA_s
                    else:
                        rth_ja, cth_ja, rth_ca, cth_ca = Rth_JA, Cth_JA, Rth_CA, Cth_CA
                    [timeTher['sw'][top.id6[j]], Tinit_T[:, j]] = calcTherRC(Tinit_T[:, j], Ta, timeLoss['sw'][top.id2[j]]['p_T'], t_ref[start:ende], rth_ja, cth_ja)
                    [timeTher['sw'][top.id8[j]], Tinit_C[:, j]] = calcTherRC(Tinit_C[:, j], Ta, timeLoss['sw'][top.id2[j]]['p_L'], t_ref[start:ende], rth_ca, cth_ca)
                    timeTher['sw'][top.id6[j]] = timeTher['sw'][top.id6[j]][:] + timeTher['sw'][top.id8[j]][:] - Ta
                elif setup['Par']['Ther']['Heatsink'] == 1 and setup['Par']['Ther']['Coupling'] == 2:
                    if setup['Top']['sourceType'] == 'DAB' and j >= 4:
                        rth_ja, cth_ja, rth_ca, cth_ca = Rth_JA_s, Cth_JA_s, Rth_CA_s, Cth_CA_s
                    else:
                        rth_ja, cth_ja, rth_ca, cth_ca = Rth_JA, Cth_JA, Rth_CA, Cth_CA
                    [timeTher['sw'][top.id6[j]], Tinit_T[:, j]] = calcTherRC(Tinit_T[:, j], Ta, timeLoss['sw'][top.id2[j]]['p_T'], t_ref[start:ende], rth_ja, cth_ja)
                    [timeTher['sw'][top.id8[j]], Tinit_C[:, j]] = calcTherRC(Tinit_C[:, j], Ta, np.mean(timeLoss['sw'][top.id2[j]]['p_L'])*np.ones(np.size(timeLoss['sw'][top.id2[j]]['p_L'])), t_ref[start:ende], rth_ca, cth_ca)
                    timeTher['sw'][top.id6[j]] = timeTher['sw'][top.id6[j]][:] + timeTher['sw'][top.id8[j]][:] - Ta
                else:
                    if setup['Top']['sourceType'] == 'DAB' and j >= 4:
                        rth_ja, cth_ja = Rth_JA_s, Cth_JA_s
                    else:
                        rth_ja, cth_ja = Rth_JA, Cth_JA
                    [timeTher['sw'][top.id6[j]], Tinit_T[:, j]] = calcTherRC(Tinit_T[:, j], Ta, timeLoss['sw'][top.id2[j]]['p_T'], t_ref[start:ende], rth_ja, cth_ja)

            # Capacitor
            timeElec['cap']['C1']['i_c'] = timeDc['i_c'][start:ende]
            timeElec['cap']['C1']['v_c'] = timeDc['v_dc'][start:ende]
            timeLoss['cap']['C1'] = calcLossCap(t_ref, timeDc['i_c'][start:ende], Tcap, para, setup)
            [timeTher['cap']['C1'], Tinit_Cap] = calcTherRC(Tinit_Cap, Ta, timeLoss['cap']['C1']['p_L'], t_ref[start:ende], Rth_JA_cap, Cth_JA_cap)

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
        if setup['Top']['sourceType'] == 'DAB' and i >= 4:
            rth_ja, cth_ja, rth_da, cth_da, rth_ca, cth_ca = Rth_JA_s, Cth_JA_s, Rth_DA_s, Cth_DA_s, Rth_CA_s, Cth_CA_s
        else:
            rth_ja, cth_ja, rth_da, cth_da, rth_ca, cth_ca = Rth_JA, Cth_JA, Rth_DA, Cth_DA, Rth_CA, Cth_CA

        [out['ther']['sw'][top.id6[i]], _] = calcTherRC(0, Ta, out['loss']['sw'][top.id2[i]]['p_T'].values, t, rth_ja, cth_ja)
        [out['ther']['sw'][top.id7[i]], _] = calcTherRC(0, Ta, out['loss']['sw'][top.id2[i]]['p_D'].values, t, rth_da, cth_da)

    # Capacitor
    [out['ther']['cap']['C1'], _] = calcTherRC(0, Ta, out['loss']['cap']['C1']['p_L'].values, t, Rth_JA_cap, Cth_JA_cap)

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
    out['ther']['cap'] = pd.DataFrame(out['ther']['cap'], columns=['C1'])

    ###################################################################################################################
    # Post-Processing
    ###################################################################################################################
    # ==============================================================================
    # Frequency domain
    # ==============================================================================
    [freqSw, freqAc, freqDc] = calcFreq(s['A'][Nsim:(K * Nsim + 1)], xs['A'][Nsim:(K * Nsim + 1)], timeAc['i_a'], timeAc['v_a'],
                                        timeAc['v_a0'], timeDc['i_dc'], timeDc['v_dc'])

    # ==============================================================================
    # Output
    # ==============================================================================
    [time, freq, _] = top.out(out['elec'], out['loss'], out['ther'], timeAc, timeDc, freqSw, freqAc, freqDc, [], [], t_ref,
                              v_ref, e_ref, s, c, xs, xsh, x, xN0, [], Nsim * (K - 1), (K * Nsim + 1), 1)

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
