#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         calcTransB4
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
from src.topo.B4.calcSSeqB4 import calcSSeqB4_CB, calcSSeqB4_FF
from src.topo.B4.calcTimeB4 import calcTimeB4
from src.general.genWaveform import genWave
from src.topo.B4.initB4 import initB4_Data, initB4
from src.topo.B4.outB4 import outB4_Trans
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
def calcTransB4(mdl, para, setupTopo, setupData, setupPara, setupExp):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("------------------------------------------")
    print("START: Transient simulation B4")
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

    # ------------------------------------------
    # IDs
    # ------------------------------------------
    id2 = ['S1', 'S2', 'S3', 'S4']
    id3 = ['A', 'A', 'B', 'B']
    id4 = [+1, +1, -1, -1]
    id5 = ['HS', 'LS', 'HS', 'LS']
    id6 = ['T1', 'T2', 'T3', 'T4']
    id7 = ['D1', 'D2', 'D3', 'D4']
    id8 = ['C1', 'C2', 'C3', 'C4']

    # ==============================================================================
    # Parameters
    # ==============================================================================
    fel = setupTopo['fel']
    fsim = setupExp['fsim']
    fs = setupPara['PWM']['fs']
    Tel = 1/fel 
    Nsim = int(np.ceil(fsim/fel))
    Npwm = int(np.ceil(fs/fel))
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
    Tj = setupData['stat']['Tj'] * np.ones((4, 1))
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
    out = initB4_Data()

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
    v_ref['A'] = +(Vdc / 2) * Mi * genWave(t_ref, fel, phiV, 0, setupTopo)
    v_ref['B'] = -(Vdc / 2) * Mi * genWave(t_ref, fel, phiV, 0, setupTopo)
    e_ref = E * genWave(t_ref, fel, phiE, 0, setupTopo) 
    
    # ==============================================================================
    # Thermal ROM
    # ==============================================================================
    [Rth_JA, Cth_JA, Rth_DA, Cth_DA, Rth_CA, Cth_CA, Rth_JA_cap, Cth_JA_cap] = initRC(para, setupPara)

    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    # ==============================================================================
    # Init
    # ==============================================================================
    Tinit_T = np.zeros((len(Rth_JA), len(id2)))
    Tinit_C = np.zeros((len(Rth_JA), len(id2)))
    Tinit_Cap = np.zeros(np.size(Rth_JA_cap))

    # ==============================================================================
    # Switching Function
    # ==============================================================================
    if setupPara['PWM']['type'] == "FF":
        [xs, xsh, s, c] = calcSSeqB4_FF(v_ref, t_ref, Mi, setupPara, setupTopo)
    elif setupPara['PWM']['type'] == "CB":
        [xs, xsh, s, c] = calcSSeqB4_CB(v_ref, t_ref, Mi, setupPara, setupTopo)
    else:
        [xs, xsh, s, c] = calcSSeqB4_CB(v_ref, t_ref, Mi, setupPara, setupTopo)

    # ==============================================================================
    # Time Domain
    # ==============================================================================
    [timeAc, timeDc] = calcTimeB4(t_ref, s, e_ref, Vdc, Mi, mdl, setupTopo, Nsim*(K-1), (K*Nsim + 1))

    # ==============================================================================
    # Electrical cycle
    # ==============================================================================
    for _ in tqdm(range(iterNel), desc='Elec-Period', position=0):
        # ------------------------------------------
        # Init
        # ------------------------------------------
        dataFel = initB4_Data()

        # ------------------------------------------
        # Electrical
        # ------------------------------------------
        timeDc['v_dc'] = calcElecCap(t_ref, timeDc['i_c'], Tcap, para, setupPara, setupTopo)

        # ------------------------------------------
        # PWM Period
        # ------------------------------------------
        for ii in range(iterNpwm):
            # Init
            [_, timeElec, timeLoss, timeTher, _, _, _, _, _] = initB4(2)
            start = int(ii*(Nsim/iterNpwm))
            if iterNpwm == 1:
                ende = int(Nsim/iterNpwm * (ii+1) + 1)
            else:
                ende = int(Nsim/iterNpwm * (ii+1) + 0)

            # Switch
            for j in range(0, len(id2)):
                timeElec['sw'][id2[j]] = calcElecSwi(Vdc, id4[j] * timeAc['i_a'][start:ende], (s[id3[j]][start:ende] == (-1) ** j), Tj[j], id5[j], para, setupPara)
                timeLoss['sw'][id2[j]] = calcLossSwi(s[id3[j]][start:ende] * (-1) ** j, timeElec['sw'][id2[j]]['i_T'], timeElec['sw'][id2[j]]['i_D'], timeElec['sw'][id2[j]]['v_T'], timeElec['sw'][id2[j]]['v_D'], Tj[j], para, setupPara)

                if setupPara['Ther']['Heatsink'] == 1 & setupPara['Ther']['Coupling'] == 1:
                    [timeTher['sw'][id6[j]], Tinit_T[:, j]] = calcTherRC(Tinit_T[:, j], Ta, timeLoss['sw'][id2[j]]['p_T'], t_ref[start:ende], Rth_JA, Cth_JA)
                    [timeTher['sw'][id8[j]], Tinit_C[:, j]] = calcTherRC(Tinit_C[:, j], Ta, timeLoss['sw'][id2[j]]['p_L'], t_ref[start:ende], Rth_CA, Cth_CA)
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
        if setupPara['Ther']['Heatsink'] == 1 & setupPara['Ther']['Coupling'] == 1:
            [out['ther']['sw'][id8[i]], _] = calcTherRC(0, Ta, out['loss']['sw'][id2[i]]['p_L'].values, t, Rth_CA, Cth_CA)
            out['ther']['sw'][id6[i]] = out['ther']['sw'][id6[i]][:] + out['ther']['sw'][id8[i]][:] - Ta
            out['ther']['sw'][id7[i]] = out['ther']['sw'][id7[i]][:] + out['ther']['sw'][id8[i]][:] - Ta
        else:
            out['ther']['sw'][id8[i]] = Ta

    # ------------------------------------------
    # Appending
    # ------------------------------------------
    # Switch
    out['ther']['sw'] = pd.DataFrame(out['ther']['sw'], columns=['T1', 'T2', 'T3', 'T4', 'D1', 'D2', 'D3', 'D4', 'C1', 'C2', 'C3', 'C4'])

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
    [time, freq] = outB4_Trans(out, timeAc, timeDc, freqSw, freqAc, freqDc, t_ref, v_ref, e_ref, s, c, xs, xsh, K, Nsim,
                               Tel, Nel)

    ###################################################################################################################
    # MSG Out
    ###################################################################################################################
    print("------------------------------------------")
    print("END: Transient simulation B4")
    print("------------------------------------------")

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return [time, freq]
