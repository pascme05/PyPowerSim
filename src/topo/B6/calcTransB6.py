#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         calcTransB6
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
from src.topo.B6.calcSSeqB6 import calcSSeqB6_CB, calcSSeqB6_FF, calcSSeqB6_SV, calcSSeqB6_OPP
from src.topo.B6.calcTimeB6 import calcTimeB6
from src.pwm.genWaveform import genWave
from src.topo.B6.initB6 import initB6_Data, initB6
from src.topo.B6.outB6 import outB6_Trans
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
def calcTransB6(mdl, para, setupTopo, setupData, setupPara, setupExp):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("------------------------------------------")
    print("START: Transient simulation B6")
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
    fs = setupPara['PWM']['fs']
    Tel = 1 / fel
    Nsim = int(np.ceil(fsim / fel))
    Npwm = int(np.ceil(fs / fel))
    K = setupData['stat']['cyc']
    Nel = int(np.ceil(setupData['trans']['tmax'] * fel))
    Mi = setupData['stat']['Mi']

    # ==============================================================================
    # Variables
    # ==============================================================================
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
    Ta = setupData['stat']['Tc']
    Tj = setupData['stat']['Tj'] * np.ones((6, 1))
    Tcap = setupData['stat']['Tj']

    # ==============================================================================
    # Update Frequency
    # ==============================================================================
    if setupExp['therFeed'] == 0:
        iterNel = 1
        iterNpwm = 1
    elif setupExp['freqPar'] == 'fel':
        iterNel = Nel
        iterNpwm = 1
    else:
        iterNel = Nel
        iterNpwm = Npwm

    # ==============================================================================
    # Outputs
    # ==============================================================================
    out = initB6_Data()

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
    for i in range(0, len(id1)):
        v_ref[id1[i]] = (Vdc / 2) * Mi * genWave(t_ref, fel, phiV, -i * 2 / 3 * np.pi, setupTopo)
        e_ref[id1[i]] = E * genWave(t_ref, fel, phiE, -i * 2 / 3 * np.pi, setupTopo)

    # ==============================================================================
    # Thermal ROM
    # ==============================================================================
    # ------------------------------------------
    # Parameters
    # ------------------------------------------
    [Rth_JA, Cth_JA, Rth_DA, Cth_DA, Rth_CA, Cth_CA, Rth_JA_cap, Cth_JA_cap] = initRC(para, setupPara)

    # ------------------------------------------
    # Init
    # ------------------------------------------
    Tinit_T = np.zeros((len(Rth_JA), len(id2)))
    Tinit_C = np.zeros((len(Rth_JA), len(id2)))
    Tinit_Cap = np.zeros(np.size(Rth_JA_cap))

    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    # ==============================================================================
    # Switching Function
    # ==============================================================================
    if setupPara['PWM']['type'] == "FF":
        [xs, xsh, s, c, x, n0] = calcSSeqB6_FF(v_ref, t_ref, Mi, setupPara, setupTopo)
    elif setupPara['PWM']['type'] == "CB":
        [xs, xsh, s, c, x, n0] = calcSSeqB6_CB(v_ref, t_ref, Mi, setupPara, setupTopo)
    elif setupPara['PWM']['type'] == "SV":
        [xs, xsh, s, c, x, n0] = calcSSeqB6_SV(v_ref, t_ref, Mi, setupPara, setupTopo)
    elif setupPara['PWM']['type'] == "OPP":
        [xs, xsh, s, c, x, n0] = calcSSeqB6_OPP(v_ref, t_ref, Mi, setupPara, setupTopo)
    else:
        [xs, xsh, s, c, x, n0] = calcSSeqB6_CB(v_ref, t_ref, Mi, setupPara, setupTopo)

    # ==============================================================================
    # Time Domain
    # ==============================================================================
    [timeAc, timeDc] = calcTimeB6(t_ref, s, e_ref, Vdc, Mi, mdl, setupTopo, Nsim * (K - 1), (K * Nsim + 1))

    # ==============================================================================
    # Electrical cycle
    # ==============================================================================
    for _ in tqdm(range(iterNel), desc='Elec-Period', position=0):
        # ------------------------------------------
        # Init
        # ------------------------------------------
        dataFel = initB6_Data()

        # ------------------------------------------
        # Electrical
        # ------------------------------------------
        timeDc['v_dc'] = calcElecCap(t_ref, timeDc['i_c'], Tcap, para, setupPara, setupTopo)

        # ------------------------------------------
        # PWM Period
        # ------------------------------------------
        for ii in range(iterNpwm):
            # Init
            [_, timeElec, timeLoss, timeTher, _, _, _, _, _] = initB6(2)
            start = int(ii * (Nsim / iterNpwm))
            if iterNpwm == 1:
                ende = int(Nsim / iterNpwm * (ii + 1) + 1)
            else:
                ende = int(Nsim / iterNpwm * (ii + 1) + 0)

            # Switch
            for j in range(0, len(id2)):
                timeElec['sw'][id2[j]] = calcElecSwi(Vdc, timeAc[id4[j]][start:ende], (s[id3[j]][start:ende] == (-1) ** j), Tj[j], id5[j], para, setupPara)
                timeLoss['sw'][id2[j]] = calcLossSwi(s[id3[j]][start:ende] * (-1) ** j, timeElec['sw'][id2[j]]['i_T'], timeElec['sw'][id2[j]]['i_D'], timeElec['sw'][id2[j]]['v_T'], timeElec['sw'][id2[j]]['v_D'], Tj[j], para, setupPara)

                if setupPara['Ther']['Heatsink'] == 1 and setupPara['Ther']['Coupling'] == 1:
                    [timeTher['sw'][id6[j]], Tinit_T[:, j]] = calcTherRC(Tinit_T[:, j], Ta, timeLoss['sw'][id2[j]]['p_T'], t_ref[start:ende], Rth_JA, Cth_JA)
                    [timeTher['sw'][id8[j]], Tinit_C[:, j]] = calcTherRC(Tinit_C[:, j], Ta, timeLoss['sw'][id2[j]]['p_L'], t_ref[start:ende], Rth_CA, Cth_CA)
                    timeTher['sw'][id6[j]] = timeTher['sw'][id6[j]][:] + timeTher['sw'][id8[j]][:] - Ta
                elif setupPara['Ther']['Heatsink'] == 1 and setupPara['Ther']['Coupling'] == 2:
                    [timeTher['sw'][id6[j]], Tinit_T[:, j]] = calcTherRC(Tinit_T[:, j], Ta, timeLoss['sw'][id2[j]]['p_T'], t_ref[start:ende], Rth_JA, Cth_JA)
                    [timeTher['sw'][id8[j]], Tinit_C[:, j]] = calcTherRC(Tinit_C[:, j], Ta, np.mean(timeLoss['sw'][id2[j]]['p_L'])*np.ones(np.size(timeLoss['sw'][id2[j]]['p_L'])), t_ref[start:ende], Rth_CA, Cth_CA)
                    timeTher['sw'][id6[j]] = timeTher['sw'][id6[j]][:] + timeTher['sw'][id8[j]][:] - Ta
                else:
                    [timeTher['sw'][id6[j]], Tinit_T[:, j]] = calcTherRC(Tinit_T[:, j], Ta, timeLoss['sw'][id2[j]]['p_T'], t_ref[start:ende], Rth_JA, Cth_JA)

            # Capacitor
            timeElec['cap']['C1']['i_c'] = timeDc['i_c'][start:ende]
            timeElec['cap']['C1']['v_c'] = timeDc['v_dc'][start:ende]
            timeLoss['cap']['C1'] = calcLossCap(t_ref, timeDc['i_c'][start:ende], Tcap, para, setupPara, setupTopo)
            [timeTher['cap']['C1'], Tinit_Cap] = calcTherRC(Tinit_Cap, Ta, timeLoss['cap']['C1']['p_L'], t_ref[start:ende], Rth_JA_cap, Cth_JA_cap)

            # Appending
            dataFel = app_fs(dataFel, timeElec, timeLoss, setupExp)

            # Parameter Update
            if setupExp['therFeed'] == 1:
                for j in range(0, len(id2)):
                    Tj[j] = timeTher['sw'][id6[j]][-1]
                Tcap = timeTher['cap']['C1'][-1]

        # ------------------------------------------
        # Appending
        # ------------------------------------------
        out = app_fel(out, dataFel['elec'], dataFel['loss'], Nel, setupExp)

    # ==============================================================================
    # Averaging
    # ==============================================================================
    out = calcAvg(out, setupExp, setupTopo, Nsim, Npwm)

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
        if setupPara['Ther']['Heatsink'] == 1 and setupPara['Ther']['Coupling'] == 1:
            [out['ther']['sw'][id8[i]], _] = calcTherRC(0, Ta, out['loss']['sw'][id2[i]]['p_L'].values, t, Rth_CA, Cth_CA)
            out['ther']['sw'][id6[i]] = out['ther']['sw'][id6[i]][:] + out['ther']['sw'][id8[i]][:] - Ta
            out['ther']['sw'][id7[i]] = out['ther']['sw'][id7[i]][:] + out['ther']['sw'][id8[i]][:] - Ta
        elif setupPara['Ther']['Heatsink'] == 1 and setupPara['Ther']['Coupling'] == 2:
            Pv_Base = np.mean((out['loss']['sw']['S1']['p_L'] + out['loss']['sw']['S2']['p_L'] + out['loss']['sw']['S3']['p_L'] +
                               out['loss']['sw']['S4']['p_L'] + out['loss']['sw']['S5']['p_L'] + out['loss']['sw']['S6']['p_L'])) / 6
            [out['ther']['sw'][id8[i]], _] = calcTherRC(0, Ta, Pv_Base*np.ones(np.size(out['loss']['sw']['S1']['p_L'])), t, Rth_CA, Cth_CA)
            out['ther']['sw'][id6[i]] = out['ther']['sw'][id6[i]][:] + out['ther']['sw'][id8[i]][:] - Ta
            out['ther']['sw'][id7[i]] = out['ther']['sw'][id7[i]][:] + out['ther']['sw'][id8[i]][:] - Ta
        else:
            out['ther']['sw'][id8[i]] = Ta

    # ------------------------------------------
    # Appending
    # ------------------------------------------
    # Switch
    out['ther']['sw'] = pd.DataFrame(out['ther']['sw'], columns=['T1', 'T2', 'T3', 'T4', 'T5', 'T6',
                                                                 'D1', 'D2', 'D3', 'D4', 'D5', 'D6',
                                                                 'C1', 'C2', 'C3', 'C4', 'C5', 'C6'])

    # Capacitor
    out['ther']['cap'] = pd.DataFrame(out['ther']['cap'], columns=['C1'])

    ###################################################################################################################
    # Post-Processing
    ###################################################################################################################
    # ==============================================================================
    # Frequency domain
    # ==============================================================================
    [freqSw, freqAc, freqDc] = calcFreq(s['A'][Nsim:(K * Nsim + 1)], xs['A'][Nsim:(K * Nsim + 1)], timeAc, timeDc)

    # ==============================================================================
    # Output
    # ==============================================================================
    [time, freq] = outB6_Trans(out, timeAc, timeDc, freqSw, freqAc, freqDc, t_ref, v_ref, e_ref, s, c, xs, xsh, x,
                               n0, K, Nsim, Tel,  Nel)

    ###################################################################################################################
    # MSG Out
    ###################################################################################################################
    print("------------------------------------------")
    print("END: Transient simulation B6")
    print("------------------------------------------")

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return [time, freq]
