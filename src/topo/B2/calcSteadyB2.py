#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         calcSteadyB2
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
from src.topo.B2.calcSSeqB2 import calcSSeqB2_CB, calcSSeqB2_FF
from src.topo.B2.calcTimeB2 import calcTimeB2
from src.general.calcFreq import calcFreq
from src.elec.calcElecSwi import calcElecSwi
from src.elec.calcLossSwi import calcLossSwi
from src.elec.calcLossCap import calcLossCap
from src.therm.calcTherRC import calcTherRC
from src.topo.B2.initB2 import initB2
from src.general.genWaveform import genWave
from src.therm.initRC import initRC
from src.elec.calcElecCap import calcElecCap
from src.topo.B2.outB2 import outB2_Steady

# ==============================================================================
# External
# ==============================================================================
import numpy as np
import math


#######################################################################################################################
# Function
#######################################################################################################################
def calcSteadyB2(mdl, para, setupTopo, setupData, setupPara, setupExp):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("------------------------------------------")
    print("START: Steady-state analysis B2 bridge")
    print("------------------------------------------")

    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Parameters
    # ==============================================================================
    fel = setupTopo['fel']
    fsim = setupExp['fsim']
    N = int(fsim/fel)
    K = setupData['stat']['cyc']
    W = setupData['stat']['W'] 
    Mi = setupData['stat']['Mi']

    # ==============================================================================
    # Variables
    # ==============================================================================
    # ------------------------------------------
    # Init 
    # ------------------------------------------
    [_, timeElec, timeLoss, timeTher, _, _, _, _, _] = initB2(W)

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
    # ------------------------------------------
    # Time
    # ------------------------------------------
    t = np.linspace(0, K / fel, K * N + 1)

    # ------------------------------------------
    # Reference
    # ------------------------------------------
    v_ref = (Vdc/2) * Mi * genWave(t, fel, phiV, 0, setupTopo)
    e_ref = E * genWave(t, fel, phiE, 0, setupTopo) 
    
    # ==============================================================================
    # Start and End
    # ==============================================================================
    start = int(N) * 2
    ende = int(K * N + 1)
    
    # ==============================================================================
    # Thermal ROM
    # ==============================================================================
    [Rth_JA, Cth_JA, Rth_DA, Cth_DA, Rth_CA, Cth_CA, Rth_JA_cap, Cth_JA_cap] = initRC(para, setupPara)

    ###################################################################################################################
    # Calculation (Stationary)
    ###################################################################################################################
    # ==============================================================================
    # Switching Function
    # ==============================================================================
    if setupPara['PWM']['type'] == "FF":
        [xs, xsh, s, c] = calcSSeqB2_FF(v_ref, t, Mi, setupPara, setupTopo)
    elif setupPara['PWM']['type'] == "CB":
        [xs, xsh, s, c] = calcSSeqB2_CB(v_ref, t, Mi, setupPara, setupTopo)
    else:
        [xs, xsh, s, c] = calcSSeqB2_CB(v_ref, t, Mi, setupPara, setupTopo)

    # ==============================================================================
    # Time Domain
    # ==============================================================================
    [timeAc, timeDc] = calcTimeB2(t, s, e_ref, Vdc, Mi, mdl, setupTopo, start, ende)

    # ==============================================================================
    # Open-Loop
    # ==============================================================================
    if setupExp['loop'] == "OL":
        # ------------------------------------------
        # Electrical
        # ------------------------------------------
        # Switches
        timeElec['sw']['S1'] = calcElecSwi(Vdc, timeAc['i_a'], (s[start:ende] == +1), Tj, 'HS', para, setupPara)
        timeElec['sw']['S2'] = calcElecSwi(Vdc, timeAc['i_a'], (s[start:ende] == -1), Tj, 'LS', para, setupPara)

        # Capacitor
        timeDc['v_dc'] = calcElecCap(t, timeDc['i_c'], Tcap, para, setupPara, setupTopo)
        timeElec['cap']['C1']['i_c'] = timeDc['i_c']
        timeElec['cap']['C1']['v_c'] = timeDc['v_dc']
        
        # ------------------------------------------
        # Losses
        # ------------------------------------------
        # Switches
        timeLoss['sw']['S1'] = calcLossSwi(s[start:ende]*(+1), timeElec['sw']['S1']['i_T'], timeElec['sw']['S1']['i_D'], timeElec['sw']['S1']['v_T'], timeElec['sw']['S1']['v_D'], Tj, para, setupPara, setupExp)
        timeLoss['sw']['S2'] = calcLossSwi(s[start:ende]*(-1), timeElec['sw']['S2']['i_T'], timeElec['sw']['S2']['i_D'], timeElec['sw']['S2']['v_T'], timeElec['sw']['S2']['v_D'], Tj, para, setupPara, setupExp)

        # Capacitor
        timeLoss['cap']['C1'] = calcLossCap(t, timeDc['i_c'], Tcap, para, setupPara, setupTopo)
        
        # ------------------------------------------
        # Thermal
        # ------------------------------------------
        [timeTher['sw']['T1'], _] = calcTherRC(0, Tc, timeLoss['sw']['S1']['p_T'], t[start:ende], Rth_JA, Cth_JA)
        [timeTher['sw']['T2'], _] = calcTherRC(0, Tc, timeLoss['sw']['S2']['p_T'], t[start:ende], Rth_JA, Cth_JA)     
        [timeTher['sw']['D1'], _] = calcTherRC(0, Tc, timeLoss['sw']['S1']['p_D'], t[start:ende], Rth_DA, Cth_DA)
        [timeTher['sw']['D2'], _] = calcTherRC(0, Tc, timeLoss['sw']['S2']['p_D'], t[start:ende], Rth_DA, Cth_DA)   
        [timeTher['cap']['C1'], _] = calcTherRC(0, Tc, timeLoss['cap']['C1']['p_L'], t[start:ende], Rth_JA_cap, Cth_JA_cap) 

    # ==============================================================================
    # Closed-Loop
    # ==============================================================================
    if setupExp['loop'] == "CL":
        # ------------------------------------------
        # Init
        # ------------------------------------------
        err = np.inf
        T_old = np.ones((3, 1))*Tj
        T_new = np.ones((3, 1))*Tj
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
            timeElec['sw']['S1'] = calcElecSwi(Vdc, timeAc['i_a'], (s[start:ende] == +1), T_old[0], 'HS', para, setupPara)
            timeElec['sw']['S2'] = calcElecSwi(Vdc, timeAc['i_a'], (s[start:ende] == -1), T_old[1], 'LS', para, setupPara)
            timeDc['v_dc'] = calcElecCap(t, timeDc['i_c'], T_old[2], para, setupPara, setupTopo)
            timeElec['cap']['C1']['i_c'] = timeDc['i_c']
            timeElec['cap']['C1']['v_c'] = timeDc['v_dc']

            # Losses
            timeLoss['sw']['S1'] = calcLossSwi(s[start:ende]*(+1), timeElec['sw']['S1']['i_T'], timeElec['sw']['S1']['i_D'], timeElec['sw']['S1']['v_T'], timeElec['sw']['S1']['v_D'], T_old[0], para, setupPara, setupExp)
            timeLoss['sw']['S2'] = calcLossSwi(s[start:ende]*(-1), timeElec['sw']['S2']['i_T'], timeElec['sw']['S2']['i_D'], timeElec['sw']['S2']['v_T'], timeElec['sw']['S2']['v_D'], T_old[1], para, setupPara, setupExp)
            timeLoss['cap']['C1'] = calcLossCap(t, timeDc['i_c'], T_old[2], para, setupPara, setupTopo)

            # Thermal
            if setupPara['Ther']['Heatsink'] == 1 & setupPara['Ther']['Coupling'] == 1:
                T_new[0] = np.mean(timeLoss['sw']['S1']['p_T']) * np.sum(Rth_JA) + np.mean(timeLoss['sw']['S1']['p_L']) * np.sum(Rth_CA) + Tc
                T_new[1] = np.mean(timeLoss['sw']['S1']['p_T']) * np.sum(Rth_JA) + np.mean(timeLoss['sw']['S1']['p_L']) * np.sum(Rth_CA) + Tc
            else:
                T_new[0] = np.mean(timeLoss['sw']['S1']['p_T']) * np.sum(Rth_JA) + Tc
                T_new[1] = np.mean(timeLoss['sw']['S2']['p_T']) * np.sum(Rth_JA) + Tc
            T_new[2] = np.mean(timeLoss['cap']['C1']['p_L']) * np.sum(Rth_JA_cap) + Tc

            # Error
            err = np.sum(abs(T_old - T_new)) / np.sum(T_old)

            # Update
            T_old = T_new[:]
            iter = iter + 1

            # Msg
            print("ITER: %d) Stationary temperature T_swi=%.2f (T_cap=%.2f) and P_swi=%.2f (Pv_cap=%.2f) with error: %.3f" % (iter, T_old[0], T_old[2], np.mean(timeLoss['sw']['S1']['p_T']), np.mean(timeLoss['cap']['C1']['p_L']), err*100))
        
        # ------------------------------------------
        # Converged
        # ------------------------------------------
        # Msg
        print("------------------------------------------")
        print("INFO: Converged after %d iterations with T_swi=%.2f (T_cap=%.2f) and error: %.3f" % (iter, T_old[0], T_old[2], err*100))

        # Thermal
        print("INFO: Calculating transient temperatures")
        [timeTher['sw']['T1'], _] = calcTherRC(0, Tc, timeLoss['sw']['S1']['p_T'], t[start:ende], Rth_JA, Cth_JA)
        [timeTher['sw']['T2'], _] = calcTherRC(0, Tc, timeLoss['sw']['S2']['p_T'], t[start:ende], Rth_JA, Cth_JA)     
        [timeTher['sw']['D1'], _] = calcTherRC(0, Tc, timeLoss['sw']['S1']['p_D'], t[start:ende], Rth_DA, Cth_DA)
        [timeTher['sw']['D2'], _] = calcTherRC(0, Tc, timeLoss['sw']['S2']['p_D'], t[start:ende], Rth_DA, Cth_DA)
        [timeTher['cap']['C1'], _] = calcTherRC(0, Tc, timeLoss['cap']['C1']['p_L'], t[start:ende], Rth_JA_cap, Cth_JA_cap)

    # ==============================================================================
    # Transient Thermal
    # ==============================================================================
    # ------------------------------------------
    # Msg
    # ------------------------------------------
    print("INFO: Calculating transient temperatures")

    # ------------------------------------------
    # Thermal
    # ------------------------------------------
    # Switches
    [timeTher['sw']['T1'], _] = calcTherRC(0, Tc, timeLoss['sw']['S1']['p_T'], t[start:ende], Rth_JA, Cth_JA)
    [timeTher['sw']['T2'], _] = calcTherRC(0, Tc, timeLoss['sw']['S2']['p_T'], t[start:ende], Rth_JA, Cth_JA)
    [timeTher['sw']['D1'], _] = calcTherRC(0, Tc, timeLoss['sw']['S1']['p_D'], t[start:ende], Rth_DA, Cth_DA)
    [timeTher['sw']['D2'], _] = calcTherRC(0, Tc, timeLoss['sw']['S2']['p_D'], t[start:ende], Rth_DA, Cth_DA)

    # Capacitor
    [timeTher['cap']['C1'], _] = calcTherRC(0, Tc, timeLoss['cap']['C1']['p_L'], t[start:ende], Rth_JA_cap, Cth_JA_cap)

    # Coupling
    if setupPara['Ther']['Heatsink'] == 1 & setupPara['Ther']['Coupling'] == 1:
        [timeTher['sw']['C1'], _] = calcTherRC(0, Tc, timeLoss['sw']['S1']['p_L'], t[start:ende], Rth_CA, Cth_CA)
        [timeTher['sw']['C2'], _] = calcTherRC(0, Tc, timeLoss['sw']['S2']['p_L'], t[start:ende], Rth_CA, Cth_CA)
        timeTher['sw']['T1'] = timeTher['sw']['T1'][:] + timeTher['sw']['C1'][:] - Tc
        timeTher['sw']['T2'] = timeTher['sw']['T2'][:] + timeTher['sw']['C2'][:] - Tc
        timeTher['sw']['D1'] = timeTher['sw']['D1'][:] + timeTher['sw']['D1'][:] - Tc
        timeTher['sw']['D2'] = timeTher['sw']['D2'][:] + timeTher['sw']['D2'][:] - Tc
    else:
        timeTher['sw']['C1'] = Tc * np.ones(np.size(s[start:ende]))
        timeTher['sw']['C2'] = Tc * np.ones(np.size(s[start:ende]))

    ###################################################################################################################
    # Post-Processing
    ###################################################################################################################
    # ==============================================================================
    # Frequency domain
    # ==============================================================================
    [freqSw, freqAc, freqDc] = calcFreq(s[start:ende], xs[start:ende], timeAc, timeDc)
    
    # ==============================================================================
    # Output
    # ==============================================================================
    [time, freq] = outB2_Steady(timeElec, timeLoss, timeTher, timeAc, timeDc, freqSw, freqAc, freqDc, t, v_ref, e_ref,
                                s, c, xs, xsh, W, ende, start)
    
    ###################################################################################################################
    # MSG Out
    ###################################################################################################################
    print("------------------------------------------")
    print("END: Steady-state analysis B2 bridge")
    print("------------------------------------------")

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return [time, freq]
