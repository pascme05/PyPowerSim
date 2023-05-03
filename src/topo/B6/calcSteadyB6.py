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
from src.topo.B6.calcFreqB6 import calcFreqB6
from src.elec.calcElecSwi import calcElecSwi
from src.elec.calcLossSwi import calcLossSwi
from src.elec.calcLossCap import calcLossCap
from src.therm.calcTherRC import calcTherRC
from src.topo.B6.initB6 import initB6
from src.general.genWaveform import genWave
from src.therm.initRC import initRC
from src.elec.calcElecCap import calcElecCap

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
    v_ref = {}
    e_ref = {}

    # ==============================================================================
    # Parameters
    # ==============================================================================
    fel = setupTopo['fel']
    fsim = setupExp['fsim']
    N = int(fsim/fel)
    K = setupData['stat']['cyc']*2
    W = setupData['stat']['W'] 
    Mi = setupData['stat']['Mi']

    # ==============================================================================
    # Variables
    # ==============================================================================
    # ------------------------------------------
    # Init 
    # ------------------------------------------
    [timeSw, timeElec, timeLoss, timeTher, freqSw, freqDc, freqAc, _, _] = initB6(W)

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
    t = np.linspace(0, K/fel, K*N+1)
    v_ref['A'] = (Vdc/2) * Mi * genWave(t, fel, phiV, 0, setupTopo)
    v_ref['B'] = (Vdc/2) * Mi * genWave(t, fel, phiV, -2/3*np.pi, setupTopo)
    v_ref['C'] = (Vdc/2) * Mi * genWave(t, fel, phiV, -4/3*np.pi, setupTopo)
    e_ref['A'] = E * genWave(t, fel, phiE, 0, setupTopo) 
    e_ref['B'] = E * genWave(t, fel, phiE, -2/3*np.pi, setupTopo) 
    e_ref['C'] = E * genWave(t, fel, phiE, -4/3*np.pi, setupTopo) 

    # ==============================================================================
    # Start and End
    # ==============================================================================
    start = int(N * (K/2))
    ende = int(K*N + 1)
    
    # ==============================================================================
    # Thermal ROM
    # ==============================================================================
    [Rth_JA, Cth_JA, Rth_DA, Cth_DA, Rth_JA_cap, Cth_JA_cap] = initRC(para, setupPara)

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

    # ==============================================================================
    # Time Domain
    # ==============================================================================
    [timeAc, timeDc] = calcTimeB6(t, s, e_ref, Vdc, Mi, mdl, setupTopo, start, ende)

    # ==============================================================================
    # Open-Loop
    # ==============================================================================
    if setupExp['loop'] == "OL":
        # ------------------------------------------
        # Electrical
        # ------------------------------------------
        timeElec['sw']['S1'] = calcElecSwi(Vdc, timeAc['i_a'], (s['A'][start:ende] ==  1), Tj, 'HS', para, setupPara)
        timeElec['sw']['S2'] = calcElecSwi(Vdc, timeAc['i_a'], (s['A'][start:ende] == -1), Tj, 'LS', para, setupPara)
        timeElec['sw']['S3'] = calcElecSwi(Vdc, timeAc['i_b'], (s['B'][start:ende] ==  1), Tj, 'HS', para, setupPara)
        timeElec['sw']['S4'] = calcElecSwi(Vdc, timeAc['i_b'], (s['B'][start:ende] == -1), Tj, 'LS', para, setupPara)
        timeElec['sw']['S5'] = calcElecSwi(Vdc, timeAc['i_c'], (s['C'][start:ende] ==  1), Tj, 'HS', para, setupPara)
        timeElec['sw']['S6'] = calcElecSwi(Vdc, timeAc['i_c'], (s['C'][start:ende] == -1), Tj, 'LS', para, setupPara)
        timeDc['v_dc'] = calcElecCap(t, timeDc['i_c'], Tcap, para, setupPara, setupTopo)
        timeElec['cap']['C1']['i_c'] = timeDc['i_c']
        timeElec['cap']['C1']['v_c'] = timeDc['v_dc']

        # ------------------------------------------
        # Losses
        # ------------------------------------------
        timeLoss['sw']['S1'] = calcLossSwi(s['A'][start:ende]*(+1), timeElec['sw']['S1']['i_T'], timeElec['sw']['S1']['i_D'], timeElec['sw']['S1']['v_T'], timeElec['sw']['S1']['v_D'], Tj, para, setupPara, setupTopo)
        timeLoss['sw']['S2'] = calcLossSwi(s['A'][start:ende]*(-1), timeElec['sw']['S2']['i_T'], timeElec['sw']['S2']['i_D'], timeElec['sw']['S2']['v_T'], timeElec['sw']['S2']['v_D'], Tj, para, setupPara, setupTopo)
        timeLoss['sw']['S3'] = calcLossSwi(s['B'][start:ende]*(+1), timeElec['sw']['S3']['i_T'], timeElec['sw']['S3']['i_D'], timeElec['sw']['S3']['v_T'], timeElec['sw']['S3']['v_D'], Tj, para, setupPara, setupTopo)
        timeLoss['sw']['S4'] = calcLossSwi(s['B'][start:ende]*(-1), timeElec['sw']['S4']['i_T'], timeElec['sw']['S4']['i_D'], timeElec['sw']['S4']['v_T'], timeElec['sw']['S4']['v_D'], Tj, para, setupPara, setupTopo)
        timeLoss['sw']['S5'] = calcLossSwi(s['C'][start:ende]*(+1), timeElec['sw']['S5']['i_T'], timeElec['sw']['S5']['i_D'], timeElec['sw']['S5']['v_T'], timeElec['sw']['S5']['v_D'], Tj, para, setupPara, setupTopo)
        timeLoss['sw']['S6'] = calcLossSwi(s['C'][start:ende]*(-1), timeElec['sw']['S6']['i_T'], timeElec['sw']['S6']['i_D'], timeElec['sw']['S6']['v_T'], timeElec['sw']['S6']['v_D'], Tj, para, setupPara, setupTopo)
        timeLoss['cap']['C1'] = calcLossCap(t, timeDc['i_c'], Tcap, para, setupPara, setupTopo)
        
        # ------------------------------------------
        # Thermal
        # ------------------------------------------
        [timeTher['sw']['T1'], _] = calcTherRC(0, Tc, timeLoss['sw']['S1']['p_T'], t[start:ende], Rth_JA, Cth_JA)
        [timeTher['sw']['T2'], _] = calcTherRC(0, Tc, timeLoss['sw']['S2']['p_T'], t[start:ende], Rth_JA, Cth_JA)
        [timeTher['sw']['T3'], _] = calcTherRC(0, Tc, timeLoss['sw']['S3']['p_T'], t[start:ende], Rth_JA, Cth_JA) 
        [timeTher['sw']['T4'], _] = calcTherRC(0, Tc, timeLoss['sw']['S4']['p_T'], t[start:ende], Rth_JA, Cth_JA)      
        [timeTher['sw']['T5'], _] = calcTherRC(0, Tc, timeLoss['sw']['S5']['p_T'], t[start:ende], Rth_JA, Cth_JA) 
        [timeTher['sw']['T6'], _] = calcTherRC(0, Tc, timeLoss['sw']['S6']['p_T'], t[start:ende], Rth_JA, Cth_JA)
        [timeTher['sw']['D1'], _] = calcTherRC(0, Tc, timeLoss['sw']['S1']['p_D'], t[start:ende], Rth_DA, Cth_DA)
        [timeTher['sw']['D2'], _] = calcTherRC(0, Tc, timeLoss['sw']['S2']['p_D'], t[start:ende], Rth_DA, Cth_DA)
        [timeTher['sw']['D3'], _] = calcTherRC(0, Tc, timeLoss['sw']['S3']['p_D'], t[start:ende], Rth_DA, Cth_DA) 
        [timeTher['sw']['D4'], _] = calcTherRC(0, Tc, timeLoss['sw']['S4']['p_D'], t[start:ende], Rth_DA, Cth_DA)    
        [timeTher['sw']['D5'], _] = calcTherRC(0, Tc, timeLoss['sw']['S5']['p_D'], t[start:ende], Rth_DA, Cth_DA) 
        [timeTher['sw']['D6'], _] = calcTherRC(0, Tc, timeLoss['sw']['S6']['p_D'], t[start:ende], Rth_DA, Cth_DA)    
        [timeTher['cap']['C1'], _] = calcTherRC(0, Tc, timeLoss['cap']['C1']['p_L'], t[start:ende], Rth_JA_cap, Cth_JA_cap) 

    # ==============================================================================
    # Open-Loop
    # ==============================================================================
    if setupExp['loop'] == "CL":
        # ------------------------------------------
        # Init
        # ------------------------------------------
        err = np.inf
        T_old = np.ones((7,1))*Tj
        T_new = np.ones((7,1))*Tj
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
            timeElec['sw']['S1'] = calcElecSwi(Vdc, timeAc['i_a'], (s['A'][start:ende] ==  1), T_old[0], 'HS', para, setupPara)
            timeElec['sw']['S2'] = calcElecSwi(Vdc, timeAc['i_a'], (s['A'][start:ende] == -1), T_old[1], 'LS', para, setupPara)
            timeElec['sw']['S3'] = calcElecSwi(Vdc, timeAc['i_b'], (s['B'][start:ende] ==  1), T_old[2], 'HS', para, setupPara)
            timeElec['sw']['S4'] = calcElecSwi(Vdc, timeAc['i_b'], (s['B'][start:ende] == -1), T_old[3], 'LS', para, setupPara)
            timeElec['sw']['S5'] = calcElecSwi(Vdc, timeAc['i_c'], (s['C'][start:ende] ==  1), T_old[4], 'HS', para, setupPara)
            timeElec['sw']['S6'] = calcElecSwi(Vdc, timeAc['i_c'], (s['C'][start:ende] == -1), T_old[5], 'LS', para, setupPara)
            timeDc['v_dc'] = calcElecCap(t, timeDc['i_c'], T_old[6], para, setupPara, setupTopo)
            timeElec['cap']['C1']['i_c'] = timeDc['i_c']
            timeElec['cap']['C1']['v_c'] = timeDc['v_dc']

            # Losses
            timeLoss['sw']['S1'] = calcLossSwi(s['A'][start:ende]*(+1), timeElec['sw']['S1']['i_T'], timeElec['sw']['S1']['i_D'], timeElec['sw']['S1']['v_T'], timeElec['sw']['S1']['v_D'], T_old[0], para, setupPara, setupTopo)
            timeLoss['sw']['S2'] = calcLossSwi(s['A'][start:ende]*(-1), timeElec['sw']['S2']['i_T'], timeElec['sw']['S2']['i_D'], timeElec['sw']['S2']['v_T'], timeElec['sw']['S2']['v_D'], T_old[1], para, setupPara, setupTopo)
            timeLoss['sw']['S3'] = calcLossSwi(s['B'][start:ende]*(+1), timeElec['sw']['S3']['i_T'], timeElec['sw']['S3']['i_D'], timeElec['sw']['S3']['v_T'], timeElec['sw']['S3']['v_D'], T_old[2], para, setupPara, setupTopo)
            timeLoss['sw']['S4'] = calcLossSwi(s['B'][start:ende]*(-1), timeElec['sw']['S4']['i_T'], timeElec['sw']['S4']['i_D'], timeElec['sw']['S4']['v_T'], timeElec['sw']['S4']['v_D'], T_old[3], para, setupPara, setupTopo)
            timeLoss['sw']['S5'] = calcLossSwi(s['C'][start:ende]*(+1), timeElec['sw']['S5']['i_T'], timeElec['sw']['S5']['i_D'], timeElec['sw']['S5']['v_T'], timeElec['sw']['S5']['v_D'], T_old[4], para, setupPara, setupTopo)
            timeLoss['sw']['S6'] = calcLossSwi(s['C'][start:ende]*(-1), timeElec['sw']['S6']['i_T'], timeElec['sw']['S6']['i_D'], timeElec['sw']['S6']['v_T'], timeElec['sw']['S6']['v_D'], T_old[5], para, setupPara, setupTopo)
            timeLoss['cap']['C1'] = calcLossCap(t, timeDc['i_c'], T_old[6], para, setupPara, setupTopo)

            # Thermal
            T_new[0] = np.mean(timeLoss['sw']['S1']['p_T']) * np.sum(Rth_JA) + Tc
            T_new[1] = np.mean(timeLoss['sw']['S2']['p_T']) * np.sum(Rth_JA) + Tc
            T_new[2] = np.mean(timeLoss['sw']['S3']['p_T']) * np.sum(Rth_JA) + Tc
            T_new[3] = np.mean(timeLoss['sw']['S4']['p_T']) * np.sum(Rth_JA) + Tc
            T_new[4] = np.mean(timeLoss['sw']['S5']['p_T']) * np.sum(Rth_JA) + Tc
            T_new[5] = np.mean(timeLoss['sw']['S6']['p_T']) * np.sum(Rth_JA) + Tc
            T_new[6] = np.mean(timeLoss['cap']['C1']['p_L']) * np.sum(Rth_JA_cap) + Tc

            # Error
            err = np.sum(abs(T_old - T_new)) / np.sum(T_old)

            # Update
            T_old = T_new[:]
            iter = iter + 1

            # Msg
            print("ITER: %d) Stationary temperatur T_swi=%.2f (T_cap=%.2f) and P_swi=%.2f (Pv_cap=%.2f) with error: %.3f" % (iter, T_old[0], T_old[6], np.mean(timeLoss['sw']['S1']['p_T']), np.mean(timeLoss['cap']['C1']['p_L']), err*100))
        
        # ------------------------------------------
        # Converged
        # ------------------------------------------
        # Msg
        print("------------------------------------------")
        print("INFO: Converged after %d iterations with T_swi=%.2f (T_cap=%.2f) and error: %.3f" % (iter, T_old[0], T_old[6], err*100))

        # Thermal
        print("INFO: Calculating transient temperatures")
        [timeTher['sw']['T1'], _] = calcTherRC(0, Tc, timeLoss['sw']['S1']['p_T'], t[start:ende], Rth_JA, Cth_JA)
        [timeTher['sw']['T2'], _] = calcTherRC(0, Tc, timeLoss['sw']['S2']['p_T'], t[start:ende], Rth_JA, Cth_JA)
        [timeTher['sw']['T3'], _] = calcTherRC(0, Tc, timeLoss['sw']['S3']['p_T'], t[start:ende], Rth_JA, Cth_JA) 
        [timeTher['sw']['T4'], _] = calcTherRC(0, Tc, timeLoss['sw']['S4']['p_T'], t[start:ende], Rth_JA, Cth_JA)      
        [timeTher['sw']['T5'], _] = calcTherRC(0, Tc, timeLoss['sw']['S5']['p_T'], t[start:ende], Rth_JA, Cth_JA) 
        [timeTher['sw']['T6'], _] = calcTherRC(0, Tc, timeLoss['sw']['S6']['p_T'], t[start:ende], Rth_JA, Cth_JA)
        [timeTher['sw']['D1'], _] = calcTherRC(0, Tc, timeLoss['sw']['S1']['p_D'], t[start:ende], Rth_DA, Cth_DA)
        [timeTher['sw']['D2'], _] = calcTherRC(0, Tc, timeLoss['sw']['S2']['p_D'], t[start:ende], Rth_DA, Cth_DA)
        [timeTher['sw']['D3'], _] = calcTherRC(0, Tc, timeLoss['sw']['S3']['p_D'], t[start:ende], Rth_DA, Cth_DA) 
        [timeTher['sw']['D4'], _] = calcTherRC(0, Tc, timeLoss['sw']['S4']['p_D'], t[start:ende], Rth_DA, Cth_DA)    
        [timeTher['sw']['D5'], _] = calcTherRC(0, Tc, timeLoss['sw']['S5']['p_D'], t[start:ende], Rth_DA, Cth_DA) 
        [timeTher['sw']['D6'], _] = calcTherRC(0, Tc, timeLoss['sw']['S6']['p_D'], t[start:ende], Rth_DA, Cth_DA) 
        [timeTher['cap']['C1'], _] = calcTherRC(0, Tc, timeLoss['cap']['C1']['p_L'], t[start:ende], Rth_JA_cap, Cth_JA_cap)
         
    ###################################################################################################################
    # Post-Processing
    ###################################################################################################################
    # ==============================================================================
    # Frequency domain
    # ==============================================================================
    [freqSw, freqAc, freqDc] = calcFreqB6(s['A'][start:ende], xs['A'][start:ende], timeAc, timeDc)
    
    # ==============================================================================
    # Output
    # ==============================================================================
    timeSw['t'] = t[0:(ende-start)]
    timeSw['v_a_ref'] = v_ref['A'][start:ende]
    timeSw['v_b_ref'] = v_ref['B'][start:ende]
    timeSw['v_c_ref'] = v_ref['C'][start:ende]
    timeSw['e_a'] = e_ref['A'][start:ende]
    timeSw['e_b'] = e_ref['B'][start:ende]
    timeSw['e_c'] = e_ref['C'][start:ende]
    timeSw['sA'] = s['A'][start:ende]
    timeSw['sB'] = s['B'][start:ende]
    timeSw['sC'] = s['C'][start:ende]
    timeSw['c'] = c[start:ende]
    timeSw['xAs'] = xs['A'][start:ende]
    timeSw['xBs'] = xs['B'][start:ende]
    timeSw['xCs'] = xs['C'][start:ende]
    timeSw['xAsh'] = xsh['A'][start:ende]
    timeSw['xBsh'] = xsh['B'][start:ende]
    timeSw['xCsh'] = xsh['C'][start:ende]
    timeSw['xA'] = x['A'][start:ende]
    timeSw['xB'] = x['B'][start:ende]
    timeSw['xC'] = x['C'][start:ende]
    timeSw['n0'] = n0[start:ende]
    
    # ==============================================================================
    # Combine
    # ==============================================================================
    time = {}
    time['Sw'] = timeSw
    time['Ac'] = timeAc
    time['Dc'] = timeDc
    time['Elec'] = timeElec
    time['Loss'] = timeLoss
    time['Ther'] = timeTher
    freq = {}
    freq['Sw'] = freqSw
    freq['Ac'] = freqAc
    freq['Dc'] = freqDc
    
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