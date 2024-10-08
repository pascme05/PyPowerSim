#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         calcClose
# Date:         14.05.2024
# Author:       Dr. Pascal A. Schirmer
# Version:      V.1.0
# Copyright:    Pascal Schirmer
#######################################################################################################################
#######################################################################################################################

#######################################################################################################################
# Function Description
#######################################################################################################################
"""
This function calculates the results for a closed loop solution for any given topology class.
Inputs:     1) top:     topology class
            2) mdl:     all models and transfer functions of the architecture
            3) para:    all parameters used in the simulation
            4) setup:   includes all simulation variables
Outputs:    1) time:    results in the time domain
            2) freq:    results in the frequency domain
"""
import copy

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
from src.elec.calcElecCap import calcElecCap

# ==============================================================================
# External
# ==============================================================================
import numpy as np
import math
from tqdm import tqdm


#######################################################################################################################
# Function
#######################################################################################################################
def calcClose(top, mdl, para, setup):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("------------------------------------------")
    print("START: Closed loop solution class", top.name)
    print("------------------------------------------")

    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Parameters
    # ==============================================================================
    fel = setup['Top']['fel']
    fsim = setup['Exp']['fsim']
    fc = setup['Par']['Cont']['fc']
    Tel = 1 / fel
    Nsim = int(np.ceil(fsim / fel))
    Ncon = int(fsim/fc)
    K = int(setup['Dat']['stat']['cyc'])
    Nel = int(np.ceil(setup['Dat']['trans']['tmax'] * fel))
    Mi = setup['Dat']['stat']['Mi']

    # ==============================================================================
    # Variables
    # ==============================================================================
    E = setup['Top']['E']
    Vdc = setup['Dat']['stat']['Vdc']
    Tj = setup['Dat']['stat']['Tj']
    phiE = math.radians(setup['Top']['phiE'])
    phiV = math.radians(setup['Dat']['stat']['phi'])

    # ==============================================================================
    # Update Frequency
    # ==============================================================================
    iterCon = int(np.ceil(setup['Dat']['trans']['tmax'] * fc))
    t_scale = np.ones(iterCon * Ncon)
    th = 0.5

    # ==============================================================================
    # Outputs
    # ==============================================================================
    [_, timeElec, timeLoss, _, _, _, _, _, _] = top.initOut()

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
    t_tot = np.linspace(0, Nel * Tel, iterCon * Ncon)

    # ------------------------------------------
    # Reference
    # ------------------------------------------
    [v_ref, e_ref, _] = top.calcRef(E, phiE, phiV, [], setup)
    [_, e_tot, i_tot] = top.calcRef(E, phiE, phiV, t_tot, setup)

    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    # ==============================================================================
    # Switching Function
    # ==============================================================================
    [xs, xsh, _, c, x, xN0] = top.calcPWM(v_ref, t_ref, Mi, setup)

    # ==============================================================================
    # Controller Init
    # ==============================================================================
    [s_i, i_act, outSw] = top.initCON()

    # ==============================================================================
    # Step Response
    # ==============================================================================
    for i in tqdm(range(iterCon), desc='Control-Periods', position=0):
        # ------------------------------------------
        # Controller
        # ------------------------------------------
        if i/iterCon < th:
            scale = 0.50
        else:
            scale = 1.00
        t_scale[i * Ncon:(i + 1) * Ncon] = scale

        # ------------------------------------------
        # Controller
        # ------------------------------------------
        [s_i, Mi_con, _] = top.calcCON(i_tot, i_act, s_i, i*Ncon-1, scale, setup)

        # ------------------------------------------
        # Append Result
        # ------------------------------------------
        [outSw, e_con, t_con] = top.appCON(s_i, outSw, i)

        # ------------------------------------------
        # Calculate Output
        # ------------------------------------------
        [tempAc, _, _] = top.calcTime(outSw, e_con, t_con, Mi_con, mdl, 0, len(t_con), [], 0, setup)
        i_act = copy.deepcopy(tempAc)

    ###################################################################################################################
    # Post-Processing
    ###################################################################################################################
    # ==============================================================================
    # Electrical Results
    # ==============================================================================
    # ------------------------------------------
    # Phase and Source
    # ------------------------------------------
    [outAc, outDc, _] = top.calcTime(outSw, e_tot, t_tot, Mi, mdl, 0, len(t_tot), [], 0, setup)
    outAc['i_ref'] = {key: value * t_scale for key, value in i_tot.items()}

    # ------------------------------------------
    # Switching Devices
    # ------------------------------------------
    for j in range(0, len(top.id2)):
        timeElec['sw'][top.id2[j]] = calcElecSwi(Vdc, top.id9[j] * outAc[top.id4[j]], (outSw[top.id3[j]] == (-1) ** j),
                                                 Tj, top.id5[j], para, setup)
        timeLoss['sw'][top.id2[j]] = calcLossSwi(outSw[top.id3[j]] * (-1) ** j,
                                                 timeElec['sw'][top.id2[j]]['i_T'], timeElec['sw'][top.id2[j]]['i_D'],
                                                 timeElec['sw'][top.id2[j]]['v_T'], timeElec['sw'][top.id2[j]]['v_D'],
                                                 Tj, para, setup)

    # ------------------------------------------
    # Capacitor
    # ------------------------------------------
    outDc['v_dc'] = calcElecCap(t_tot, outDc['i_c'], Tj, para, setup)
    timeLoss['cap']['C1'] = calcLossCap(t_tot, outDc['i_c'], Tj, para, setup)
    timeElec['cap']['C1']['v_c'] = outDc['v_dc']
    timeElec['cap']['C1']['i_c'] = outDc['i_c']

    # ==============================================================================
    # Frequency domain
    # ==============================================================================
    [freqSw, freqAc, freqDc] = calcFreq(outSw['A'][-(K * Nsim + 1):-1], xs['A'][-(K * Nsim + 1):-1],
                                        outAc['i_a'][-(K * Nsim + 1):-1], outAc['v_a'][-(K * Nsim + 1):-1],
                                        outAc['v_a0'][-(K * Nsim + 1):-1], outDc['i_dc'][-(K * Nsim + 1):-1],
                                        outDc['v_dc'][-(K * Nsim + 1):-1])

    # ==============================================================================
    # Output
    # ==============================================================================
    [time, freq, _] = top.out(timeElec, timeLoss, [], outAc, outDc, freqSw, freqAc, freqDc, [], [], t_ref,
                              v_ref, e_ref, outSw, c, xs, xsh, x, xN0, [], 0, Nsim + 1, 1)

    ###################################################################################################################
    # MSG Out
    ###################################################################################################################
    print("------------------------------------------")
    print("END: Closed loop class", top.name)
    print("------------------------------------------")

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return [time, freq]
