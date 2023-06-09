#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         calcTransB4
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
from src.topo.B4.calcSSeqB4 import calcSSeqB4_CB, calcSSeqB4_FF
from src.topo.B4.calcTimeB4 import calcTimeB4
from src.general.genWaveform import genWave
from src.topo.B4.initB4 import initB4_Data, initB4
from src.topo.B4.calcFreqB4 import calcFreqB4
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
    v_ref = {}

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
    Tj = [setupData['stat']['Tj'], setupData['stat']['Tj'], setupData['stat']['Tj'], setupData['stat']['Tj']]
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
    out = initB4_Data()

    ###################################################################################################################
    # Pre-Processing
    ###################################################################################################################
    # ==============================================================================
    # Generate Reference Waveform
    # ==============================================================================
    t_ref = np.linspace(0, K/fel, K*Nsim+1)
    v_ref['A'] = +(Vdc/2) * Mi * genWave(t_ref, fel, phiV, 0, setupTopo)
    v_ref['B'] = -(Vdc/2) * Mi * genWave(t_ref, fel, phiV, 0, setupTopo)
    e_ref = E * genWave(t_ref, fel, phiE, 0, setupTopo) 
    
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
    Tinit_C1 = np.zeros(np.size(Rth_JA_cap))

    # ==============================================================================
    # Switching Function
    # ==============================================================================
    if setupPara['PWM']['type'] == "FF":
        [xs, xsh, s, c] = calcSSeqB4_FF(v_ref, t_ref, Mi, setupPara, setupTopo)
    elif setupPara['PWM']['type'] == "CB":
        [xs, xsh, s, c] = calcSSeqB4_CB(v_ref, t_ref, Mi, setupPara, setupTopo)

    # ==============================================================================
    # Time Domain
    # ==============================================================================
    [timeAc, timeDc] = calcTimeB4(t_ref, s, e_ref, Vdc, Mi, mdl, setupTopo, Nsim*(K-1), (K*Nsim + 1))

    # ==============================================================================
    # Elelctrical cycle
    # ==============================================================================
    for i in trange(iterNel, desc='Elec-Period'):
        # ------------------------------------------
        # Init
        # ------------------------------------------
        dataFel = initB4_Data()
        [_, timeElec, timeLoss, timeTher, _, _, _, _, _] = initB4(2)

        # ------------------------------------------
        # Electrical
        # ------------------------------------------
        timeDc['v_dc'] = calcElecCap(t_ref, timeDc['i_c'], Tcap, para, setupPara, setupTopo)

        # ------------------------------------------
        # PWM Period
        # ------------------------------------------
        for ii in trange(iterNpwm, desc='PWM-Period', leave=False):
            # Init
            [_, timeElec, timeLoss, timeTher, _, _, _, _, _] = initB4(2)
            start = int(ii*(Nsim/iterNpwm))
            if iterNpwm == 1:
                ende = int(Nsim/iterNpwm * (ii+1) + 1)
            else:
                ende = int(Nsim/iterNpwm * (ii+1) + 0)

            # Electrical
            timeElec['sw']['S1'] = calcElecSwi(Vdc,  timeAc['i_a'][start:ende], (s['A'][start:ende] ==  1), Tj[0], 'HS', para, setupPara)
            timeElec['sw']['S2'] = calcElecSwi(Vdc,  timeAc['i_a'][start:ende], (s['A'][start:ende] == -1), Tj[1], 'LS', para, setupPara)
            timeElec['sw']['S3'] = calcElecSwi(Vdc, -timeAc['i_a'][start:ende], (s['B'][start:ende] ==  1), Tj[2], 'HS', para, setupPara)
            timeElec['sw']['S4'] = calcElecSwi(Vdc, -timeAc['i_a'][start:ende], (s['B'][start:ende] == -1), Tj[3], 'LS', para, setupPara)
            timeElec['cap']['C1']['i_c'] = timeDc['i_c'][start:ende]
            timeElec['cap']['C1']['v_c'] = timeDc['v_dc'][start:ende]

            # Losses
            timeLoss['sw']['S1'] = calcLossSwi(s['A'][start:ende]*(+1), timeElec['sw']['S1']['i_T'], timeElec['sw']['S1']['i_D'], timeElec['sw']['S1']['v_T'], timeElec['sw']['S1']['v_D'], Tj[0], para, setupPara, setupExp)
            timeLoss['sw']['S2'] = calcLossSwi(s['A'][start:ende]*(-1), timeElec['sw']['S2']['i_T'], timeElec['sw']['S2']['i_D'], timeElec['sw']['S2']['v_T'], timeElec['sw']['S2']['v_D'], Tj[1], para, setupPara, setupExp)
            timeLoss['sw']['S3'] = calcLossSwi(s['B'][start:ende]*(+1), timeElec['sw']['S3']['i_T'], timeElec['sw']['S3']['i_D'], timeElec['sw']['S3']['v_T'], timeElec['sw']['S3']['v_D'], Tj[2], para, setupPara, setupExp)
            timeLoss['sw']['S4'] = calcLossSwi(s['B'][start:ende]*(-1), timeElec['sw']['S4']['i_T'], timeElec['sw']['S4']['i_D'], timeElec['sw']['S4']['v_T'], timeElec['sw']['S4']['v_D'], Tj[3], para, setupPara, setupExp)
            timeLoss['cap']['C1'] = calcLossCap(t_ref, timeDc['i_c'][start:ende], Tcap, para, setupPara, setupTopo)

            # Thermal
            [timeTher['sw']['T1'], Tinit_T1] = calcTherRC(Tinit_T1, Tc, timeLoss['sw']['S1']['p_T'], t_ref[start:ende], Rth_JA, Cth_JA)
            [timeTher['sw']['T2'], Tinit_T2] = calcTherRC(Tinit_T2, Tc, timeLoss['sw']['S2']['p_T'], t_ref[start:ende], Rth_JA, Cth_JA) 
            [timeTher['sw']['T3'], Tinit_T3] = calcTherRC(Tinit_T3, Tc, timeLoss['sw']['S3']['p_T'], t_ref[start:ende], Rth_JA, Cth_JA) 
            [timeTher['sw']['T4'], Tinit_T4] = calcTherRC(Tinit_T4, Tc, timeLoss['sw']['S4']['p_T'], t_ref[start:ende], Rth_JA, Cth_JA)      
            [timeTher['cap']['C1'], Tinit_C1] = calcTherRC(Tinit_C1, Tc, timeLoss['cap']['C1']['p_L'], t_ref[start:ende], Rth_JA_cap, Cth_JA_cap) 

            # Apending
            dataFel = app_fs(dataFel, timeElec, timeLoss, setupExp)

            # Parameter Update
            if setupExp['loop'] == 'CL':
                Tj[0] = timeTher['sw']['T1'][-1]
                Tj[1] = timeTher['sw']['T2'][-1]
                Tj[2] = timeTher['sw']['T3'][-1]
                Tj[3] = timeTher['sw']['T4'][-1]
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
    [out['ther']['sw']['D1'], _] = calcTherRC(0, Tc, out['loss']['sw']['S1']['p_D'].values, t, Rth_DA, Cth_DA)
    [out['ther']['sw']['D2'], _] = calcTherRC(0, Tc, out['loss']['sw']['S2']['p_D'].values, t, Rth_DA, Cth_DA)
    [out['ther']['sw']['D3'], _] = calcTherRC(0, Tc, out['loss']['sw']['S3']['p_D'].values, t, Rth_DA, Cth_DA)
    [out['ther']['sw']['D4'], _] = calcTherRC(0, Tc, out['loss']['sw']['S4']['p_D'].values, t, Rth_DA, Cth_DA)
    [out['ther']['cap']['C1'], _] = calcTherRC(0, Tc, out['loss']['cap']['C1']['p_L'].values, t, Rth_JA_cap, Cth_JA_cap)
    out['ther']['sw'] = pd.DataFrame(out['ther']['sw'], columns = ['T1','T2','T3','T4','D1','D2','D3','D4'])  
    out['ther']['cap'] = pd.DataFrame(out['ther']['cap'], columns = ['C1']) 


    ###################################################################################################################
    # Post-Processing
    ###################################################################################################################
    # ==============================================================================
    # Variables
    # ==============================================================================
    [timeSw, _, _, _, freqSw, freqDc, freqAc, _, _] = initB4(5)

    # ==============================================================================
    # Frequency domain
    # ==============================================================================
    [freqSw, freqAc, freqDc] = calcFreqB4(s['A'][Nsim:(K*Nsim + 1)], xs['A'][Nsim:(K*Nsim + 1)], timeAc, timeDc)
    
    # ==============================================================================
    # Output
    # ==============================================================================
    timeSw['t'] = t_ref[0:((K*Nsim + 1)-Nsim)]
    timeSw['v_a_ref'] = v_ref['A'][Nsim:(K*Nsim + 1)]
    timeSw['v_b_ref'] = v_ref['B'][Nsim:(K*Nsim + 1)]
    timeSw['e'] = e_ref[Nsim:(K*Nsim + 1)]
    timeSw['sA'] = s['A'][Nsim:(K*Nsim + 1)]
    timeSw['sB'] = s['B'][Nsim:(K*Nsim + 1)]
    timeSw['cA'] = c['A'][Nsim:(K*Nsim + 1)]
    timeSw['cB'] = c['B'][Nsim:(K*Nsim + 1)]
    timeSw['xAs'] = xs['A'][Nsim:(K*Nsim + 1)]
    timeSw['xBs'] = xs['B'][Nsim:(K*Nsim + 1)]
    timeSw['xAsh'] = xsh['A'][Nsim:(K*Nsim + 1)]
    timeSw['xBsh'] = xsh['B'][Nsim:(K*Nsim + 1)]
    
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
    print("END: Transient simulation B4")
    print("------------------------------------------")

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return [time, freq]