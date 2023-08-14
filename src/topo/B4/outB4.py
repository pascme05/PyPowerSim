#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         outB4
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
from src.topo.B4.initB4 import initB4

# ==============================================================================
# External
# ==============================================================================
import numpy as np


#######################################################################################################################
# Function
#######################################################################################################################
def outB4_Trans(out, timeAc, timeDc, freqSw, freqAc, freqDc, t_ref, v_ref, e_ref, s, c, xs, xsh, K, Nsim, Tel, Nel):
    ###################################################################################################################
    # Init Values
    ###################################################################################################################
    # ==============================================================================
    # Variables
    # ==============================================================================
    time = {}
    freq = {}

    # ==============================================================================
    # Struct
    # ==============================================================================
    [timeSw, _, _, _, _, _, _, _, _] = initB4(4)

    ###################################################################################################################
    # Output Values
    ###################################################################################################################
    # ==============================================================================
    # General
    # ==============================================================================
    timeSw['t'] = t_ref[0:((K * Nsim + 1) - Nsim)]
    timeSw['v_a_ref'] = v_ref['A'][Nsim:(K * Nsim + 1)]
    timeSw['v_b_ref'] = v_ref['B'][Nsim:(K * Nsim + 1)]
    timeSw['e'] = e_ref[Nsim:(K * Nsim + 1)]
    timeSw['sA'] = s['A'][Nsim:(K * Nsim + 1)]
    timeSw['sB'] = s['B'][Nsim:(K * Nsim + 1)]
    timeSw['cA'] = c['A'][Nsim:(K * Nsim + 1)]
    timeSw['cB'] = c['B'][Nsim:(K * Nsim + 1)]
    timeSw['xAs'] = xs['A'][Nsim:(K * Nsim + 1)]
    timeSw['xBs'] = xs['B'][Nsim:(K * Nsim + 1)]
    timeSw['xAsh'] = xsh['A'][Nsim:(K * Nsim + 1)]
    timeSw['xBsh'] = xsh['B'][Nsim:(K * Nsim + 1)]

    # ==============================================================================
    # Combine
    # ==============================================================================
    # ------------------------------------------
    # Time
    # ------------------------------------------
    time['Sw'] = timeSw
    time['Ac'] = timeAc
    time['Dc'] = timeDc
    time['t'] = np.linspace(0, Tel * Nel, int(len(out['loss']['sw']['S1']['p_T'])))
    time['Elec'] = out['elec']
    time['Loss'] = out['loss']
    time['Ther'] = out['ther']

    # ------------------------------------------
    # Frequency
    # ------------------------------------------
    freq['Sw'] = freqSw
    freq['Ac'] = freqAc
    freq['Dc'] = freqDc

    ###################################################################################################################
    # Outputs
    ###################################################################################################################
    return [time, freq]


#######################################################################################################################
# Function
#######################################################################################################################
def outB4_Steady(timeElec, timeLoss, timeTher, timeAc, timeDc, freqSw, freqAc, freqDc, t, v_ref, e_ref, s, c, xs, xsh, W, ende, start):
    ###################################################################################################################
    # Init Values
    ###################################################################################################################
    # ==============================================================================
    # Variables
    # ==============================================================================
    time = {}
    freq = {}

    # ==============================================================================
    # Struct
    # ==============================================================================
    [timeSw, _, _, _, _, _, _, _, _] = initB4(W)

    ###################################################################################################################
    # Output Values
    ###################################################################################################################
    # ==============================================================================
    # General
    # ==============================================================================
    timeSw['t'] = t[0:(ende - start)]
    timeSw['v_a_ref'] = v_ref['A'][start:ende]
    timeSw['v_b_ref'] = v_ref['B'][start:ende]
    timeSw['e'] = e_ref[start:ende]
    timeSw['sA'] = s['A'][start:ende]
    timeSw['sB'] = s['B'][start:ende]
    timeSw['cA'] = c['A'][start:ende]
    timeSw['cB'] = c['B'][start:ende]
    timeSw['xAs'] = xs['A'][start:ende]
    timeSw['xBs'] = xs['B'][start:ende]
    timeSw['xAsh'] = xsh['A'][start:ende]
    timeSw['xBsh'] = xsh['B'][start:ende]

    # ==============================================================================
    # Combine
    # ==============================================================================
    # ------------------------------------------
    # Time
    # ------------------------------------------
    time['Sw'] = timeSw
    time['Ac'] = timeAc
    time['Dc'] = timeDc
    time['Elec'] = timeElec
    time['Loss'] = timeLoss
    time['Ther'] = timeTher

    # ------------------------------------------
    # Frequency
    # ------------------------------------------
    freq['Sw'] = freqSw
    freq['Ac'] = freqAc
    freq['Dc'] = freqDc

    ###################################################################################################################
    # Outputs
    ###################################################################################################################
    return [time, freq]


#######################################################################################################################
# Function
#######################################################################################################################
def outB4_Sweep(distAc, distDc, timeAc, timeDc, freqSw, freqAc, freqDc, t, v_ref, e_ref, s, c, xs, xsh, W, M_i, ende, start):
    ###################################################################################################################
    # Init Values
    ###################################################################################################################
    # ==============================================================================
    # Variables
    # ==============================================================================
    time = {}
    freq = {}
    sweep = {}

    # ==============================================================================
    # Struct
    # ==============================================================================
    [timeSw, _, _, _, _, _, _, _, _] = initB4(W)

    ###################################################################################################################
    # Output Values
    ###################################################################################################################
    # ==============================================================================
    # General
    # ==============================================================================
    timeSw['t'] = t[0:(ende - start)]
    timeSw['v_a_ref'] = v_ref['A'][start:ende]
    timeSw['v_b_ref'] = v_ref['B'][start:ende]
    timeSw['e'] = e_ref[start:ende]
    timeSw['sA'] = s['A'][start:ende]
    timeSw['sB'] = s['B'][start:ende]
    timeSw['cA'] = c['A'][start:ende]
    timeSw['cB'] = c['B'][start:ende]
    timeSw['xAs'] = xs['A'][start:ende]
    timeSw['xBs'] = xs['B'][start:ende]
    timeSw['xAsh'] = xsh['A'][start:ende]
    timeSw['xBsh'] = xsh['B'][start:ende]

    # ==============================================================================
    # Combine
    # ==============================================================================
    # ------------------------------------------
    # Time
    # ------------------------------------------
    time['Sw'] = timeSw
    time['Ac'] = timeAc
    time['Dc'] = timeDc

    # ------------------------------------------
    # Frequency
    # ------------------------------------------
    freq['Sw'] = freqSw
    freq['Ac'] = freqAc
    freq['Dc'] = freqDc

    # ------------------------------------------
    # Sweep
    # ------------------------------------------
    sweep['Ac'] = distAc
    sweep['Dc'] = distDc
    sweep['Mi'] = M_i

    ###################################################################################################################
    # Outputs
    ###################################################################################################################
    return [time, freq, sweep]