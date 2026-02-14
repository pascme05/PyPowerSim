#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         calcSweep_DCDC
# Date:         04.05.2024
# Author:       Dr. Pascal A. Schirmer
# Version:      V.1.0
# Copyright:    Pascal Schirmer
#######################################################################################################################
#######################################################################################################################

#######################################################################################################################
# Function Description
#######################################################################################################################
"""
This function calculates the results for a parameter sweep for any given topology class.
Inputs:     1) top:     topology class
            2) mdl:     all models and transfer functions of the architecture
            3) para:    all parameters used in the simulation
            4) setup:   includes all simulation variables
Outputs:    1) time:    results in the time domain
            2) freq:    results in the frequency domain
            3) dist:    results in the distortion domain
"""

#######################################################################################################################
# Import libs
#######################################################################################################################
# ==============================================================================
# Internal
# ==============================================================================
from src.general.calcSpec import calcFreq as calcFreqSpec
from src.general.calcFreq import calcFreq as calcFreqStd
from src.general.calcDistNum import calcDistNum
from src.general.helpFnc import calcDistSignals

# ==============================================================================
# External
# ==============================================================================
import numpy as np
import math
from tqdm import tqdm


#######################################################################################################################
# Function
#######################################################################################################################
def calcSweep_DCDC(top, mdl, para, setup):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("------------------------------------------")
    print("START: Sweeping class", top.name)
    print("------------------------------------------")

    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Parameters
    # ==============================================================================
    fel = setup['Top']['fel']
    fsim = setup['Exp']['fsim']
    N = int(fsim / fel)
    K = int(setup['Dat']['stat']['cyc'])
    W = int(setup['Dat']['stat']['W'])
    M_i = setup['Dat']['stat']['Mi']
    E = setup['Top']['E']
    Vdc = setup['Dat']['stat']['Vdc']
    phiE = math.radians(setup['Top']['phiE'])
    phiV = math.radians(setup['Dat']['stat']['phi'])
    L = setup['Top']['L']
    R = setup['Top']['R']
    Z = np.sqrt(R ** 2 + (2 * np.pi * fel * L) ** 2)

    # ==============================================================================
    # Variables
    # ==============================================================================
    [_, _, _, _, _, _, _, distAc, distDc] = top.initOut()

    ###################################################################################################################
    # Pre-Processing
    ###################################################################################################################
    # ==============================================================================
    # Generate Reference Waveform
    # ==============================================================================
    t_ref = np.linspace(0, K / fel, K * N + 1)
    [v_ref, e_ref, _] = top.calcRef(E, phiE, phiV, [], setup)

    # ==============================================================================
    # Sweep Parameter
    # ==============================================================================
    if setup['Top']['sourceType'] == 'DAB':
        phi_base = setup['Dat']['stat'].get('PhiDAB', np.pi / 2)
        M_i = top.Mi
        # If PhiDAB is still in degrees, convert
        if abs(phi_base) > 2 * np.pi:
            phi_base = math.radians(phi_base)
        if phi_base == 0:
            phi_base = np.pi / 2
        phi_sign = 1 if phi_base >= 0 else -1
        phi_span = max(abs(phi_base), 1e-3)
        phi_i = np.linspace(1e-3, phi_span - 1e-3, W) * phi_sign
    else:
        print("ERROR: Topology does not exist.")

    # ==============================================================================
    # Start and End
    # ==============================================================================
    start = int(N) * 2
    ende = int(K * N + 1)

    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    # ==============================================================================
    # Stationary
    # ==============================================================================
    # ------------------------------------------
    # Switching Function
    # ------------------------------------------
    if setup['Top']['sourceType'] == 'DAB':
        # Use base phase shift for stationary point
        phi_base = setup['Dat']['stat'].get('PhiDAB', np.pi / 2)
        if abs(phi_base) > 2 * np.pi:
            phi_base = math.radians(phi_base)
        top.phi = phi_base
    else:
        print("ERROR: Topology does not exist.")
    [xs, xsh, s, c, x, xN0] = top.calcPWM(v_ref, t_ref, M_i, setup)

    # ------------------------------------------
    # Time Domain
    # ------------------------------------------
    if setup['Top']['sourceType'] == 'DAB':
        [timeAc, timeDc, _] = top.calcTime(s, None, None, e_ref, t_ref, M_i, mdl, start, ende, [], para, setup)
    else:
        print("ERROR: Topology does not exist.")

    # ==============================================================================
    # Sweeping DAB
    # ==============================================================================
    if setup['Top']['sourceType'] == 'DAB':
        for i in tqdm(range(len(phi_i)), desc='Sweep'):
            # ------------------------------------------
            # Switching
            # ------------------------------------------
            top.phi = phi_i[i]
            [_, _, s_i, _, _, _] = top.calcPWM(v_ref, t_ref, M_i, setup)

            # ------------------------------------------
            # Time
            # ------------------------------------------
            [tempAc, tempDc, _] = top.calcTime(s_i, None, None, e_ref, t_ref, M_i, mdl, start, ende, [], para, setup)

            # ------------------------------------------
            # Distortion
            # ------------------------------------------
            # Primary
            Vdc_sec = np.mean(tempDc['v_dc_sec'])
            dist_map = calcDistSignals(
                t_ref[start:ende],
                {
                    'i_ac_pri': tempAc['i_ac_pri'],
                    'v_ac_pri': tempAc['v_ac_pri'],
                    'i_ac_sec': tempAc['i_ac_sec'],
                    'v_ac_sec': tempAc['v_ac_sec'],
                    'i_dc_pri': tempDc['i_dc_pri'],
                    'v_dc_pri': tempDc['v_dc_pri'],
                    'i_dc_sec': tempDc['i_dc_sec'],
                    'v_dc_sec': tempDc['v_dc_sec']
                },
                f1_map={
                    'i_ac_pri': fel, 'v_ac_pri': fel,
                    'i_ac_sec': fel, 'v_ac_sec': fel,
                    'i_dc_pri': 0.0, 'v_dc_pri': 0.0,
                    'i_dc_sec': 0.0, 'v_dc_sec': 0.0
                },
                dc_offset_map={'v_dc_pri': Vdc, 'v_dc_sec': Vdc_sec},
                default_f1=fel,
                default_dc=0.0
            )
            dist_i_pri = dist_map['i_ac_pri']
            dist_v_pri = dist_map['v_ac_pri']
            dist_i_sec = dist_map['i_ac_sec']
            dist_v_sec = dist_map['v_ac_sec']
            dist_i_dc_pri = dist_map['i_dc_pri']
            dist_v_dc_pri = dist_map['v_dc_pri']
            dist_i_dc_sec = dist_map['i_dc_sec']
            dist_v_dc_sec = dist_map['v_dc_sec']
            [anaAcPri, anaDcPri] = top.calcDist(tempAc['i_ac_pri'], tempAc['v_ac_pri'], M_i, L, Z, setup)
            [anaAcSec, anaDcSec] = top.calcDist(tempAc['i_ac_sec'], tempAc['v_ac_sec'], M_i, L, Z, setup)

            distAc['Pri']['num']['V_a_eff'][i] = dist_v_pri['eff']
            distAc['Pri']['num']['V_a_v1_eff'][i] = dist_v_pri['v1_eff']
            distAc['Pri']['num']['V_a_thd'][i] = dist_v_pri['thd']
            distAc['Pri']['num']['I_a_eff'][i] = dist_i_pri['eff']
            distAc['Pri']['num']['I_a_v1_eff'][i] = dist_i_pri['v1_eff']
            distAc['Pri']['num']['I_a_thd'][i] = dist_i_pri['thd']

            distDc['Pri']['num']['V_dc_eff'][i] = dist_v_dc_pri['eff']
            distDc['Pri']['num']['V_dc_v1_eff'][i] = dist_v_dc_pri['v1_eff']
            distDc['Pri']['num']['V_dc_thd'][i] = dist_v_dc_pri['thd']
            distDc['Pri']['num']['I_dc_eff'][i] = dist_i_dc_pri['eff']
            distDc['Pri']['num']['I_dc_v1_eff'][i] = dist_i_dc_pri['v1_eff']
            distDc['Pri']['num']['I_dc_thd'][i] = dist_i_dc_pri['thd']

            distAc['Sec']['num']['V_a_eff'][i] = dist_v_sec['eff']
            distAc['Sec']['num']['V_a_v1_eff'][i] = dist_v_sec['v1_eff']
            distAc['Sec']['num']['V_a_thd'][i] = dist_v_sec['thd']
            distAc['Sec']['num']['I_a_eff'][i] = dist_i_sec['eff']
            distAc['Sec']['num']['I_a_v1_eff'][i] = dist_i_sec['v1_eff']
            distAc['Sec']['num']['I_a_thd'][i] = dist_i_sec['thd']

            distDc['Sec']['num']['V_dc_eff'][i] = dist_v_dc_sec['eff']
            distDc['Sec']['num']['V_dc_v1_eff'][i] = dist_v_dc_sec['v1_eff']
            distDc['Sec']['num']['V_dc_thd'][i] = dist_v_dc_sec['thd']
            distDc['Sec']['num']['I_dc_eff'][i] = dist_i_dc_sec['eff']
            distDc['Sec']['num']['I_dc_v1_eff'][i] = dist_i_dc_sec['v1_eff']
            distDc['Sec']['num']['I_dc_thd'][i] = dist_i_dc_sec['thd']

            for c1 in anaAcPri:
                distAc['Pri']['ana'][c1][i] = anaAcPri[c1]
            for c1 in anaDcPri:
                distDc['Pri']['ana'][c1][i] = anaDcPri[c1]
            for c1 in anaAcSec:
                distAc['Sec']['ana'][c1][i] = anaAcSec[c1]
            for c1 in anaDcSec:
                distDc['Sec']['ana'][c1][i] = anaDcSec[c1]
    else:
        print("ERROR: Topology does not exist.")

    ###################################################################################################################
    # Post-Processing
    ###################################################################################################################
    # ==============================================================================
    # Frequency domain
    # ==============================================================================
    if setup['Top']['sourceType'] == 'DAB':
        dictSw = {'Sa': s['A'][start:ende], 'Xas': xs['A'][start:ende], 'Sb': s['C'][start:ende]}
        [freqSw, freqAc, freqDc] = calcFreqSpec(dictSw, timeAc, timeDc)
    else:
        print("ERROR: Topology does not exist.")

    # ==============================================================================
    # Output
    # ==============================================================================
    [time, freq, sweep] = top.out([], [], [], timeAc, timeDc, freqSw, freqAc, freqDc, distAc, distDc, t_ref, v_ref,
                                  e_ref, s, c, xs, xsh, x, xN0, M_i, start, ende, 1)

    ###################################################################################################################
    # MSG Out
    ###################################################################################################################
    print("------------------------------------------")
    print("END: Sweeping class", top.name)
    print("------------------------------------------")

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return [time, freq, sweep]
