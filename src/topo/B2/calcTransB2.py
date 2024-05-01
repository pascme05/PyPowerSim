#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         calcTransB2
# Date:         01.05.2024
# Author:       Dr. Pascal A. Schirmer
# Version:      V.1.0
# Copyright:    Pascal Schirmer
#######################################################################################################################
#######################################################################################################################

#######################################################################################################################
# Function Description
#######################################################################################################################
"""
This function calculates the results for a transient solution of the B2 half-bridge circuit.
Inputs:     1) mdl:     all models and transfer functions of the architecture
            2) para:    all parameters used in the simulation
            3) setup:   includes all simulation variables
Outputs:    1) time:    results in the time domain
            2) freq:    results in the frequency domain
"""

#######################################################################################################################
# Import libs
#######################################################################################################################
# ==============================================================================
# Internal
# ==============================================================================
from src.topo.B2.calcSSeqB2 import calcSSeqB2_CB, calcSSeqB2_FF, calcSSeqB2_OPP
from src.topo.B2.calcTimeB2 import calcTimeB2
from src.pwm.genWaveform import genWave
from src.topo.B2.initB2 import initB2_Data
from src.general.calcFreq import calcFreq
from src.elec.calcElecSwi import calcElecSwi
from src.elec.calcLossSwi import calcLossSwi
from src.therm.calcTherRC import calcTherRC
from src.elec.calcLossCap import calcLossCap
from src.general.calcAvg import calcAvg
from src.topo.B2.initB2 import initB2
from src.therm.initRC import initRC
from src.elec.calcElecCap import calcElecCap
from src.general.append import app_fs, app_fel
from src.topo.B2.outB2 import outB2_Trans

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
def calcTransB2(mdl, para, setup):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("------------------------------------------")
    print("START: Transient simulation B2")
    print("------------------------------------------")

    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Init
    # ==============================================================================
    id2 = ['S1', 'S2']
    id5 = ['HS', 'LS']
    id6 = ['T1', 'T2']
    id7 = ['D1', 'D2']
    id8 = ['C1', 'C2']

    # ==============================================================================
    # Parameters
    # ==============================================================================
    fel = setup['Top']['fel']
    fsim = setup['Exp']['fsim']
    fs = setup['Par']['PWM']['fs']
    Tel = 1/fel 
    Nsim = int(np.ceil(fsim/fel))
    Npwm = int(np.ceil(fs/fel))
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
    Tj = [setup['Dat']['stat']['Tj'], setup['Dat']['stat']['Tj']]
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
    out = initB2_Data()

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
    v_ref = (Vdc/2) * Mi * genWave(t_ref, fel, phiV, setup)
    e_ref = E * genWave(t_ref, fel, phiE, setup)
    
    # ==============================================================================
    # Thermal ROM
    # ==============================================================================
    [Rth_JA, Cth_JA, Rth_DA, Cth_DA, Rth_CA, Cth_CA, Rth_JA_cap, Cth_JA_cap] = initRC(para, setup)

    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    # ==============================================================================
    # Init
    # ==============================================================================
    Tinit_T = np.zeros((len(Rth_JA), 2))
    Tinit_C = np.zeros((len(Rth_JA), 2))
    Tinit_Cap = np.zeros(np.size(Rth_JA_cap))

    # ==============================================================================
    # Switching Function
    # ==============================================================================
    if setup['Par']['PWM']['type'] == "FF":
        [xs, xsh, s, c] = calcSSeqB2_FF(v_ref, t_ref, Mi, setup)
    elif setup['Par']['PWM']['type'] == "CB":
        [xs, xsh, s, c] = calcSSeqB2_CB(v_ref, t_ref, Mi, setup)
    elif setup['Par']['PWM']['type'] == "OPP":
        [xs, xsh, s, c] = calcSSeqB2_OPP(v_ref, t_ref, Mi, setup)
    else:
        [xs, xsh, s, c] = calcSSeqB2_CB(v_ref, t_ref, Mi, setup)

    # ==============================================================================
    # Time Domain
    # ==============================================================================
    [timeAc, timeDc] = calcTimeB2(t_ref, s, e_ref, Vdc, Mi, mdl, setup, Nsim*(K-1), (K*Nsim + 1))

    # ==============================================================================
    # Electrical cycle
    # ==============================================================================
    for _ in tqdm(range(iterNel), desc='Elec-Period', position=0):
        # ------------------------------------------
        # Init
        # ------------------------------------------
        dataFel = initB2_Data()

        # ------------------------------------------
        # Electrical
        # ------------------------------------------
        timeDc['v_dc'] = calcElecCap(t_ref, timeDc['i_c'], Tcap, para, setup)

        # ------------------------------------------
        # PWM Period
        # ------------------------------------------
        for ii in range(iterNpwm):
            # Init
            [_, timeElec, timeLoss, timeTher, _, _, _, _, _] = initB2(2)
            start = int(ii*(Nsim/iterNpwm))
            if iterNpwm == 1:
                ende = int(Nsim/iterNpwm * (ii+1) + 1)
            else:
                ende = int(Nsim/iterNpwm * (ii+1) + 0)

            # Switch
            for j in range(0, len(id2)):
                timeElec['sw'][id2[j]] = calcElecSwi(Vdc, timeAc['i_a'][start:ende], (s[start:ende] == (-1) ** j), Tj[j], id5[j], para, setup)
                timeLoss['sw'][id2[j]] = calcLossSwi(s[start:ende] * (-1) ** j, timeElec['sw'][id2[j]]['i_T'], timeElec['sw'][id2[j]]['i_D'], timeElec['sw'][id2[j]]['v_T'], timeElec['sw'][id2[j]]['v_D'], Tj[j], para, setup)

                if setup['Par']['Ther']['Heatsink'] == 1 and setup['Par']['Ther']['Coupling'] == 1:
                    [timeTher['sw'][id6[j]], Tinit_T[:, j]] = calcTherRC(Tinit_T[:, j], Ta, timeLoss['sw'][id2[j]]['p_T'], t_ref[start:ende], Rth_JA, Cth_JA)
                    [timeTher['sw'][id8[j]], Tinit_C[:, j]] = calcTherRC(Tinit_C[:, j], Ta, timeLoss['sw'][id2[j]]['p_L'], t_ref[start:ende], Rth_CA, Cth_CA)
                    timeTher['sw'][id6[j]] = timeTher['sw'][id6[j]][:] + timeTher['sw'][id8[j]][:] - Ta
                elif setup['Par']['Ther']['Heatsink'] == 1 and setup['Par']['Ther']['Coupling'] == 2:
                    [timeTher['sw'][id6[j]], Tinit_T[:, j]] = calcTherRC(Tinit_T[:, j], Ta, timeLoss['sw'][id2[j]]['p_T'], t_ref[start:ende], Rth_JA, Cth_JA)
                    [timeTher['sw'][id8[j]], Tinit_C[:, j]] = calcTherRC(Tinit_C[:, j], Ta, np.mean(timeLoss['sw'][id2[j]]['p_L']) * np.ones(np.size(timeLoss['sw'][id2[j]]['p_L'])), t_ref[start:ende], Rth_CA, Cth_CA)
                    timeTher['sw'][id6[j]] = timeTher['sw'][id6[j]][:] + timeTher['sw'][id8[j]][:] - Ta
                else:
                    [timeTher['sw'][id6[j]], Tinit_T[:, j]] = calcTherRC(Tinit_T[:, j], Ta, timeLoss['sw'][id2[j]]['p_T'], t_ref[start:ende], Rth_JA, Cth_JA)

            # Capacitor
            timeElec['cap']['C1']['i_c'] = timeDc['i_c'][start:ende]
            timeElec['cap']['C1']['v_c'] = timeDc['v_dc'][start:ende]
            timeLoss['cap']['C1'] = calcLossCap(t_ref, timeDc['i_c'][start:ende], Tcap, para, setup)
            [timeTher['cap']['C1'], Tinit_Cap] = calcTherRC(Tinit_Cap, Ta, timeLoss['cap']['C1']['p_L'], t_ref[start:ende], Rth_JA_cap, Cth_JA_cap)

            # Appending
            dataFel = app_fs(dataFel, timeElec, timeLoss, setup)

            # Parameter Update
            if setup['Exp']['therFeed'] == 1:
                for j in range(0, len(id2)):
                    Tj[j] = timeTher['sw'][id6[j]][-1]
                Tcap = timeTher['cap']['C1'][-1]

        # ------------------------------------------
        # Appending
        # ------------------------------------------
        out = app_fel(out, dataFel['elec'], dataFel['loss'], Nel, setup)

    # ==============================================================================
    # Averaging
    # ==============================================================================
    out = calcAvg(out, setup, Nsim, Npwm)

    # ==============================================================================
    # Thermal
    # ==============================================================================
    # ------------------------------------------
    # Init
    # ------------------------------------------
    t = np.linspace(0, Tel * Nel, len(out['loss']['sw']['S1']['p_T']))

    # ------------------------------------------
    # Calc
    # ------------------------------------------
    # Switches
    for i in range(0, len(id2)):
        [out['ther']['sw'][id6[i]], _] = calcTherRC(0, Ta, out['loss']['sw'][id2[i]]['p_T'].values, t, Rth_JA, Cth_JA)
        [out['ther']['sw'][id7[i]], _] = calcTherRC(0, Ta, out['loss']['sw'][id2[i]]['p_D'].values, t, Rth_DA, Cth_DA)

    # Capacitor
    [out['ther']['cap']['C1'], _] = calcTherRC(0, Ta, out['loss']['cap']['C1']['p_L'].values, t, Rth_JA_cap, Cth_JA_cap)

    # Coupling
    for i in range(0, len(id2)):
        if setup['Par']['Ther']['Heatsink'] == 1 and setup['Par']['Ther']['Coupling'] == 1:
            [out['ther']['sw'][id8[i]], _] = calcTherRC(0, Ta, out['loss']['sw'][id2[i]]['p_L'].values, t, Rth_CA, Cth_CA)
            out['ther']['sw'][id6[i]] = out['ther']['sw'][id6[i]][:] + out['ther']['sw'][id8[i]][:] - Ta
            out['ther']['sw'][id7[i]] = out['ther']['sw'][id7[i]][:] + out['ther']['sw'][id8[i]][:] - Ta
        elif setup['Par']['Ther']['Heatsink'] == 1 and setup['Par']['Ther']['Coupling'] == 2:
            Pv_Base = np.mean((out['loss']['sw']['S1']['p_L'] + out['loss']['sw']['S2']['p_L'])) / 2
            [out['ther']['sw'][id8[i]], _] = calcTherRC(0, Ta, Pv_Base*np.ones(np.size(out['loss']['sw']['S1']['p_L'])), t, Rth_CA, Cth_CA)
            out['ther']['sw'][id6[i]] = out['ther']['sw'][id6[i]][:] + out['ther']['sw'][id8[i]][:] - Ta
            out['ther']['sw'][id7[i]] = out['ther']['sw'][id7[i]][:] + out['ther']['sw'][id8[i]][:] - Ta
        else:
            out['ther']['sw'][id8[i]] = Ta

    # ------------------------------------------
    # Appending
    # ------------------------------------------
    # Switch
    out['ther']['sw'] = pd.DataFrame(out['ther']['sw'], columns=['T1', 'T2', 'D1', 'D2', 'C1', 'C2'])

    # Capacitor
    out['ther']['cap'] = pd.DataFrame(out['ther']['cap'], columns=['C1'])

    ###################################################################################################################
    # Post-Processing
    ###################################################################################################################
    # ==============================================================================
    # Frequency domain
    # ==============================================================================
    [freqSw, freqAc, freqDc] = calcFreq(s[Nsim:(K*Nsim + 1)], xs[Nsim:(K*Nsim + 1)], timeAc, timeDc)
    
    # ==============================================================================
    # Output
    # ==============================================================================
    [time, freq] = outB2_Trans(out, timeAc, timeDc, freqSw, freqAc, freqDc, t_ref, v_ref, e_ref, s, c, xs, xsh, K, Nsim,
                               Tel, Nel)

    ###################################################################################################################
    # MSG Out
    ###################################################################################################################
    print("------------------------------------------")
    print("END: Transient simulation B2")
    print("------------------------------------------")

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return [time, freq]
