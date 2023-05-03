#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         calcTransB6
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
from src.general.genWaveform import genWave
from src.topo.B6.initB6 import initB6_Data, initB6
from src.topo.B6.calcFreqB6 import calcFreqB6
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
from tqdm.auto import trange

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
    v_ref = {}
    e_ref = {}

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
    Tc = setupData['stat']['Tc']
    Tj = setupData['stat']['Tj'] * np.ones((6,1))
    Tcap = setupData['stat']['Tj']
    
    # ==============================================================================
    # Update Frequency
    # ==============================================================================
    if setupExp['loop'] == 'OL':
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
    t_ref = np.linspace(0, K/fel, K*Nsim+1)
    v_ref['A'] = (Vdc/2) * Mi * genWave(t_ref, fel, phiV, 0, setupTopo)
    v_ref['B'] = (Vdc/2) * Mi * genWave(t_ref, fel, phiV, -2/3*np.pi, setupTopo)
    v_ref['C'] = (Vdc/2) * Mi * genWave(t_ref, fel, phiV, -4/3*np.pi, setupTopo)
    e_ref['A'] = E * genWave(t_ref, fel, phiE, 0, setupTopo) 
    e_ref['B'] = E * genWave(t_ref, fel, phiE, -2/3*np.pi, setupTopo) 
    e_ref['C'] = E * genWave(t_ref, fel, phiE, -4/3*np.pi, setupTopo) 
    
    # ==============================================================================
    # Thermal ROM
    # ==============================================================================
    [Rth_JA, Cth_JA, Rth_DA, Cth_DA, Rth_JA_cap, Cth_JA_cap] = initRC(para, setupPara)

    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    # ==============================================================================
    # Init
    # ==============================================================================
    Tinit_T1 = np.zeros(np.size(Rth_JA))
    Tinit_T2 = np.zeros(np.size(Rth_JA))
    Tinit_T3 = np.zeros(np.size(Rth_JA))
    Tinit_T4 = np.zeros(np.size(Rth_JA))
    Tinit_T5 = np.zeros(np.size(Rth_JA))
    Tinit_T6 = np.zeros(np.size(Rth_JA))
    Tinit_C1 = np.zeros(np.size(Rth_JA_cap))

    # ==============================================================================
    # Switching Function
    # ==============================================================================
    if setupPara['PWM']['type'] == "FF":
        [xs, xsh, s, c, x, n0] = calcSSeqB6_FF(v_ref, t_ref, Mi, setupPara, setupTopo)
    elif setupPara['PWM']['type'] == "CB":
        [xs, xsh, s, c, x, n0] = calcSSeqB6_CB(v_ref, t_ref, Mi, setupPara, setupTopo)
    elif setupPara['PWM']['type'] == "SV":
        [xs, xsh, s, c, x, n0] = calcSSeqB6_SV(v_ref, t_ref, Mi, setupPara, setupTopo)

    # ==============================================================================
    # Time Domain
    # ==============================================================================
    [timeAc, timeDc] = calcTimeB6(t_ref, s, e_ref, Vdc, Mi, mdl, setupTopo, Nsim*(K-1), (K*Nsim + 1))

    # ==============================================================================
    # Elelctrical cycle
    # ==============================================================================
    for i in trange(iterNel, desc='Elec-Period'):
        # ------------------------------------------
        # Init
        # ------------------------------------------
        dataFel = initB6_Data()
        [_, timeElec, timeLoss, timeTher, _, _, _, _, _] = initB6(2)

        # ------------------------------------------
        # Electrical
        # ------------------------------------------
        timeDc['v_dc'] = calcElecCap(t_ref, timeDc['i_c'], Tcap, para, setupPara, setupTopo)

        # ------------------------------------------
        # PWM Period
        # ------------------------------------------
        for ii in trange(iterNpwm, desc='PWM-Period', leave=False):
            # Init
            [_, timeElec, timeLoss, timeTher, _, _, _, _, _] = initB6(2)
            start = int(ii*(Nsim/iterNpwm))
            if iterNpwm == 1:
                ende = int(Nsim/iterNpwm * (ii+1) + 1)
            else:
                ende = int(Nsim/iterNpwm * (ii+1) + 0)

            # Electrical
            timeElec['sw']['S1'] = calcElecSwi(Vdc, timeAc['i_a'][start:ende], (s['A'][start:ende] ==  1), Tj[0], 'HS', para, setupPara)
            timeElec['sw']['S2'] = calcElecSwi(Vdc, timeAc['i_a'][start:ende], (s['A'][start:ende] == -1), Tj[1], 'LS', para, setupPara)
            timeElec['sw']['S3'] = calcElecSwi(Vdc, timeAc['i_b'][start:ende], (s['B'][start:ende] ==  1), Tj[2], 'HS', para, setupPara)
            timeElec['sw']['S4'] = calcElecSwi(Vdc, timeAc['i_b'][start:ende], (s['B'][start:ende] == -1), Tj[3], 'LS', para, setupPara)
            timeElec['sw']['S5'] = calcElecSwi(Vdc, timeAc['i_c'][start:ende], (s['C'][start:ende] ==  1), Tj[4], 'HS', para, setupPara)
            timeElec['sw']['S6'] = calcElecSwi(Vdc, timeAc['i_c'][start:ende], (s['C'][start:ende] == -1), Tj[5], 'LS', para, setupPara)
            timeElec['cap']['C1']['i_c'] = timeDc['i_c'][start:ende]
            timeElec['cap']['C1']['v_c'] = timeDc['v_dc'][start:ende]

            # Losses
            timeLoss['sw']['S1'] = calcLossSwi(s['A'][start:ende]*(+1), timeElec['sw']['S1']['i_T'], timeElec['sw']['S1']['i_D'], timeElec['sw']['S1']['v_T'], timeElec['sw']['S1']['v_D'], Tj[0], para, setupPara, setupTopo)
            timeLoss['sw']['S2'] = calcLossSwi(s['A'][start:ende]*(-1), timeElec['sw']['S2']['i_T'], timeElec['sw']['S2']['i_D'], timeElec['sw']['S2']['v_T'], timeElec['sw']['S2']['v_D'], Tj[1], para, setupPara, setupTopo)
            timeLoss['sw']['S3'] = calcLossSwi(s['B'][start:ende]*(+1), timeElec['sw']['S3']['i_T'], timeElec['sw']['S3']['i_D'], timeElec['sw']['S3']['v_T'], timeElec['sw']['S3']['v_D'], Tj[2], para, setupPara, setupTopo)
            timeLoss['sw']['S4'] = calcLossSwi(s['B'][start:ende]*(-1), timeElec['sw']['S4']['i_T'], timeElec['sw']['S4']['i_D'], timeElec['sw']['S4']['v_T'], timeElec['sw']['S4']['v_D'], Tj[3], para, setupPara, setupTopo)
            timeLoss['sw']['S5'] = calcLossSwi(s['C'][start:ende]*(+1), timeElec['sw']['S5']['i_T'], timeElec['sw']['S5']['i_D'], timeElec['sw']['S5']['v_T'], timeElec['sw']['S5']['v_D'], Tj[4], para, setupPara, setupTopo)
            timeLoss['sw']['S6'] = calcLossSwi(s['C'][start:ende]*(-1), timeElec['sw']['S6']['i_T'], timeElec['sw']['S6']['i_D'], timeElec['sw']['S6']['v_T'], timeElec['sw']['S6']['v_D'], Tj[5], para, setupPara, setupTopo)
            timeLoss['cap']['C1'] = calcLossCap(t_ref, timeDc['i_c'][start:ende], Tcap, para, setupPara, setupTopo)

            # Thermal
            [timeTher['sw']['T1'], Tinit_T1] = calcTherRC(Tinit_T1, Tc, timeLoss['sw']['S1']['p_T'], t_ref[start:ende], Rth_JA, Cth_JA)
            [timeTher['sw']['T2'], Tinit_T2] = calcTherRC(Tinit_T2, Tc, timeLoss['sw']['S2']['p_T'], t_ref[start:ende], Rth_JA, Cth_JA) 
            [timeTher['sw']['T3'], Tinit_T3] = calcTherRC(Tinit_T3, Tc, timeLoss['sw']['S3']['p_T'], t_ref[start:ende], Rth_JA, Cth_JA) 
            [timeTher['sw']['T4'], Tinit_T4] = calcTherRC(Tinit_T4, Tc, timeLoss['sw']['S4']['p_T'], t_ref[start:ende], Rth_JA, Cth_JA)      
            [timeTher['sw']['T5'], Tinit_T5] = calcTherRC(Tinit_T5, Tc, timeLoss['sw']['S5']['p_T'], t_ref[start:ende], Rth_JA, Cth_JA) 
            [timeTher['sw']['T6'], Tinit_T6] = calcTherRC(Tinit_T6, Tc, timeLoss['sw']['S6']['p_T'], t_ref[start:ende], Rth_JA, Cth_JA) 
            [timeTher['cap']['C1'], Tinit_C1] = calcTherRC(Tinit_C1, Tc, timeLoss['cap']['C1']['p_L'], t_ref[start:ende], Rth_JA_cap, Cth_JA_cap) 

            # Apending
            dataFel = app_fs(dataFel, timeElec, timeLoss, setupExp)

            # Parameter Update
            if setupExp['loop'] == 'CL':
                Tj[0] = timeTher['sw']['T1'][-1]
                Tj[1] = timeTher['sw']['T2'][-1]
                Tj[2] = timeTher['sw']['T3'][-1]
                Tj[3] = timeTher['sw']['T4'][-1]
                Tj[4] = timeTher['sw']['T5'][-1]
                Tj[5] = timeTher['sw']['T6'][-1]
                Tcap  = timeTher['cap']['C1'][-1]

        # ------------------------------------------
        # Apending
        # ------------------------------------------
        out = app_fel(out, dataFel['elec'], dataFel['loss'], Nel, setupExp)

    # ==============================================================================
    # Averaging
    # ==============================================================================
    out = calcAvg(out, setupExp, setupTopo, Nsim, Npwm)

    # ==============================================================================
    # Thermal
    # ==============================================================================
    t = np.linspace(0, Tel*Nel, len(out['loss']['sw']['S1']['p_T']))
    [out['ther']['sw']['T1'], _] = calcTherRC(0, Tc, out['loss']['sw']['S1']['p_T'].values, t, Rth_JA, Cth_JA)
    [out['ther']['sw']['T2'], _] = calcTherRC(0, Tc, out['loss']['sw']['S2']['p_T'].values, t, Rth_JA, Cth_JA)     
    [out['ther']['sw']['T3'], _] = calcTherRC(0, Tc, out['loss']['sw']['S3']['p_T'].values, t, Rth_JA, Cth_JA)
    [out['ther']['sw']['T4'], _] = calcTherRC(0, Tc, out['loss']['sw']['S4']['p_T'].values, t, Rth_JA, Cth_JA) 
    [out['ther']['sw']['T5'], _] = calcTherRC(0, Tc, out['loss']['sw']['S5']['p_T'].values, t, Rth_JA, Cth_JA)
    [out['ther']['sw']['T6'], _] = calcTherRC(0, Tc, out['loss']['sw']['S6']['p_T'].values, t, Rth_JA, Cth_JA) 
    [out['ther']['sw']['D1'], _] = calcTherRC(0, Tc, out['loss']['sw']['S1']['p_D'].values, t, Rth_DA, Cth_DA)
    [out['ther']['sw']['D2'], _] = calcTherRC(0, Tc, out['loss']['sw']['S2']['p_D'].values, t, Rth_DA, Cth_DA)
    [out['ther']['sw']['D3'], _] = calcTherRC(0, Tc, out['loss']['sw']['S3']['p_D'].values, t, Rth_DA, Cth_DA)
    [out['ther']['sw']['D4'], _] = calcTherRC(0, Tc, out['loss']['sw']['S4']['p_D'].values, t, Rth_DA, Cth_DA)
    [out['ther']['sw']['D5'], _] = calcTherRC(0, Tc, out['loss']['sw']['S5']['p_D'].values, t, Rth_DA, Cth_DA)
    [out['ther']['sw']['D6'], _] = calcTherRC(0, Tc, out['loss']['sw']['S6']['p_D'].values, t, Rth_DA, Cth_DA)
    [out['ther']['cap']['C1'], _] = calcTherRC(0, Tc, out['loss']['cap']['C1']['p_L'].values, t, Rth_JA_cap, Cth_JA_cap)
    out['ther']['sw'] = pd.DataFrame(out['ther']['sw'], columns = ['T1','T2','T3','T4','T5','T6','D1','D2','D3','D4','D5','D6'])  
    out['ther']['cap'] = pd.DataFrame(out['ther']['cap'], columns = ['C1']) 

    ###################################################################################################################
    # Post-Processing
    ###################################################################################################################
    # ==============================================================================
    # Variables
    # ==============================================================================
    [timeSw, _, _, _, freqSw, freqDc, freqAc, _, _] = initB6(5)

    # ==============================================================================
    # Frequency domain
    # ==============================================================================
    [freqSw, freqAc, freqDc] = calcFreqB6(s['A'][Nsim:(K*Nsim + 1)], xs['A'][Nsim:(K*Nsim + 1)], timeAc, timeDc)
    
    # ==============================================================================
    # Output
    # ==============================================================================
    timeSw['t'] = t_ref[0:((K*Nsim + 1)-Nsim)]
    timeSw['v_a_ref'] = v_ref['A'][Nsim:(K*Nsim + 1)]
    timeSw['v_b_ref'] = v_ref['B'][Nsim:(K*Nsim + 1)]
    timeSw['v_c_ref'] = v_ref['C'][Nsim:(K*Nsim + 1)]
    timeSw['e_a'] = e_ref['A'][Nsim:(K*Nsim + 1)]
    timeSw['e_b'] = e_ref['B'][Nsim:(K*Nsim + 1)]
    timeSw['e_c'] = e_ref['C'][Nsim:(K*Nsim + 1)]
    timeSw['sA'] = s['A'][Nsim:(K*Nsim + 1)]
    timeSw['sB'] = s['B'][Nsim:(K*Nsim + 1)]
    timeSw['sC'] = s['C'][Nsim:(K*Nsim + 1)]
    timeSw['c'] = c[Nsim:(K*Nsim + 1)]
    timeSw['xAs'] = xs['A'][Nsim:(K*Nsim + 1)]
    timeSw['xBs'] = xs['B'][Nsim:(K*Nsim + 1)]
    timeSw['xCs'] = xs['C'][Nsim:(K*Nsim + 1)]
    timeSw['xAsh'] = xsh['A'][Nsim:(K*Nsim + 1)]
    timeSw['xBsh'] = xsh['B'][Nsim:(K*Nsim + 1)]
    timeSw['xCsh'] = xsh['C'][Nsim:(K*Nsim + 1)]
    timeSw['xA'] = x['A'][Nsim:(K*Nsim + 1)]
    timeSw['xB'] = x['B'][Nsim:(K*Nsim + 1)]
    timeSw['xC'] = x['C'][Nsim:(K*Nsim + 1)]
    timeSw['n0'] = n0[Nsim:(K*Nsim + 1)]
    
    # ==============================================================================
    # Combine
    # ==============================================================================
    # ------------------------------------------
    # Time
    # ------------------------------------------
    time = {}
    time['Sw'] = timeSw
    time['Ac'] = timeAc
    time['Dc'] = timeDc
    time['t'] = np.linspace(0, Tel*Nel, int(len(out['loss']['sw']['S1']['p_T'])))
    time['Elec'] = out['elec']
    time['Loss'] = out['loss']
    time['Ther'] = out['ther']
    
    # ------------------------------------------
    # Frequency
    # ------------------------------------------
    freq = {}
    freq['Sw'] = freqSw
    freq['Ac'] = freqAc
    freq['Dc'] = freqDc

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