#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         plotResults_DCDC
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
from src.general.helpFnc import rms

# ==============================================================================
# External
# ==============================================================================
import numpy as np


#######################################################################################################################
# Function
#######################################################################################################################
def plotResults_DCDC(time, setup):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("------------------------------------------")
    print("START: Printing Results")
    print("------------------------------------------")

    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Labels
    # ==============================================================================
    def _sort_ids(ids):
        def _key(x):
            try:
                return int(x[1:])
            except:
                return 0
        return sorted(ids, key=_key)

    if 'Elec' in time and 'sw' in time['Elec']:
        id1 = _sort_ids(list(time['Elec']['sw'].keys()))
    else:
        if setup['Top']['sourceType'] == 'DAB':
            id1 = ['S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8']
        else:
            print("WARN: Invalid topology assuming DAB")
            id1 = ['S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8']

    id2 = ['T' + s[1:] for s in id1]
    id3 = ['D' + s[1:] for s in id1]
    lab = id1 + ['C1']

    # ==============================================================================
    # Variables
    # ==============================================================================
    pt_sec = time['Dc']['i_dc_sec'] * time['Dc']['v_dc_sec']
    # pt_pri = time['Dc']['i_dc_pri'] * time['Dc']['v_dc_pri']

    # ==============================================================================
    # Init
    # ==============================================================================
    I_ALL = np.zeros((len(id1) + 2 + 1, 3))  # +2 for caps, +1 for transformer
    V_ALL = np.zeros((len(id1) + 2 + 1, 3))
    P_ALL = np.zeros((len(id1) + 2 + 1, 4))
    T_ALL = np.zeros((len(id1) + 2 + 1, 4))
    E_ALL = np.zeros((len(id1) + 2 + 1, 1))

    ###################################################################################################################
    # Pre-Processing
    ###################################################################################################################
    # ==============================================================================
    # Component Level
    # ==============================================================================
    # ------------------------------------------
    # Switches
    # ------------------------------------------
    for i in range(0, len(id1)):
        # Currents
        I_ALL[i, 0] = np.max(time['Elec']['sw'][id1[i]]['i_T'] + time['Elec']['sw'][id1[i]]['i_D'])
        I_ALL[i, 1] = np.mean(time['Elec']['sw'][id1[i]]['i_T'] + time['Elec']['sw'][id1[i]]['i_D'])
        I_ALL[i, 2] = rms(time['Elec']['sw'][id1[i]]['i_T'] + time['Elec']['sw'][id1[i]]['i_D'])

        # Voltages
        V_ALL[i, 0] = np.max(time['Elec']['sw'][id1[i]]['v_T'])
        V_ALL[i, 1] = np.mean(time['Elec']['sw'][id1[i]]['v_T'])
        V_ALL[i, 2] = rms(time['Elec']['sw'][id1[i]]['v_T'])

        # Losses
        P_ALL[i, 0] = np.mean(time['Loss']['sw'][id1[i]]['p_T_s'])
        P_ALL[i, 1] = np.mean(time['Loss']['sw'][id1[i]]['p_T_c'])
        P_ALL[i, 2] = np.mean(time['Loss']['sw'][id1[i]]['p_D_s'])
        P_ALL[i, 3] = np.mean(time['Loss']['sw'][id1[i]]['p_D_c'])

        # Thermal
        T_ALL[i, 0] = np.max(time['Ther']['sw'][id2[i]])
        T_ALL[i, 1] = np.mean(time['Ther']['sw'][id2[i]])
        T_ALL[i, 2] = np.max(time['Ther']['sw'][id3[i]])
        T_ALL[i, 3] = np.mean(time['Ther']['sw'][id3[i]])

        # Efficiency
        E_ALL[i, 0] = np.abs(rms(pt_sec) - np.mean(time['Loss']['sw'][id1[i]]['p_T'])) / rms(pt_sec)

    # ------------------------------------------
    # Capacitor
    # ------------------------------------------
    # Currents
    I_ALL[len(id1), 0] = np.max(time['Elec']['cap']['C1']['i_c'])
    I_ALL[len(id1), 1] = np.mean(time['Elec']['cap']['C1']['i_c'])
    I_ALL[len(id1), 2] = rms(time['Elec']['cap']['C1']['i_c'])
    I_ALL[len(id1) + 1, 0] = np.max(time['Elec']['cap']['C2']['i_c'])
    I_ALL[len(id1) + 1, 1] = np.mean(time['Elec']['cap']['C2']['i_c'])
    I_ALL[len(id1) + 1, 2] = rms(time['Elec']['cap']['C2']['i_c'])

    # Voltages
    V_ALL[len(id1), 0] = np.max(time['Elec']['cap']['C1']['v_c'])
    V_ALL[len(id1), 1] = np.mean(time['Elec']['cap']['C1']['v_c'])
    V_ALL[len(id1), 2] = rms(time['Elec']['cap']['C1']['v_c'])
    V_ALL[len(id1) + 1, 0] = np.max(time['Elec']['cap']['C2']['v_c'])
    V_ALL[len(id1) + 1, 1] = np.mean(time['Elec']['cap']['C2']['v_c'])
    V_ALL[len(id1) + 1, 2] = rms(time['Elec']['cap']['C2']['v_c'])

    # Losses
    P_ALL[len(id1), 0] = np.mean(time['Loss']['cap']['C1']['p_L'])
    P_ALL[len(id1) + 1, 0] = np.mean(time['Loss']['cap']['C2']['p_L'])

    # Thermal
    T_ALL[len(id1), 0] = np.max(time['Ther']['cap']['C1'])
    T_ALL[len(id1), 1] = np.mean(time['Ther']['cap']['C1'])
    T_ALL[len(id1) + 1, 0] = np.max(time['Ther']['cap']['C2'])
    T_ALL[len(id1) + 1, 1] = np.mean(time['Ther']['cap']['C2'])

    # ------------------------------------------
    # Transformer
    # ------------------------------------------
    # Currents
    I_ALL[len(id1) + 2, 0] = np.max(time['Elec']['tra']['i_p'])
    I_ALL[len(id1) + 2, 1] = np.mean(time['Elec']['tra']['i_p'])
    I_ALL[len(id1) + 2, 2] = rms(time['Elec']['tra']['i_p'])

    # Voltages
    V_ALL[len(id1) + 2, 0] = np.max(time['Elec']['tra']['v_p'])
    V_ALL[len(id1) + 2, 1] = np.mean(time['Elec']['tra']['v_p'])
    V_ALL[len(id1) + 2, 2] = rms(time['Elec']['tra']['v_p'])

    # Losses
    P_ALL[len(id1) + 2, 0] = np.mean(time['Loss']['tra']['p_PC'])
    P_ALL[len(id1) + 2, 1] = np.mean(time['Loss']['tra']['p_SC'])
    P_ALL[len(id1) + 2, 2] = np.mean(time['Loss']['tra']['p_CC'])
    P_ALL[len(id1) + 2, 3] = np.mean(time['Loss']['tra']['p_L'])

    # Thermal
    T_ALL[len(id1) + 2, 0] = np.max(time['Ther']['tra']['PC'])
    T_ALL[len(id1) + 2, 1] = np.max(time['Ther']['tra']['SC'])
    T_ALL[len(id1) + 2, 2] = np.max(time['Ther']['tra']['CC'])
    T_ALL[len(id1) + 2, 3] = np.mean(time['Ther']['tra']['PC'])

    # ==============================================================================
    # Converter Level
    # ==============================================================================
    # ------------------------------------------
    # Currents and Voltages
    # ------------------------------------------
    # I_pri_rms = rms(time['Ac']['i_ac_pri'])
    I_sec_rms = rms(time['Ac']['i_ac_sec'])
    # V_pri_rms = rms(time['Ac']['v_ac_pri'])
    V_sec_rms = rms(time['Ac']['v_ac_sec'])
    I_dc_in_rms = rms(time['Dc']['i_dc_in'])
    I_dc_out_rms = rms(time['Dc']['i_dc_out'])
    V_dc_in_rms = rms(time['Dc']['v_dc_pri'])
    V_dc_out_rms = rms(time['Dc']['v_dc_sec'])

    # ------------------------------------------
    # Power
    # ------------------------------------------
    Pin = I_dc_in_rms * V_dc_in_rms
    Po = I_dc_out_rms * V_dc_out_rms

    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    def _print_overview():
        print("CONVERTER OVERVIEW")
        print('|          |                     CURRENT AND VOLTAGE                   |             TOTAL LOSSES              |     TOTAL POWER AND EFF.    |')
        print('| item ID  |  I_MAX  |  I_AVG  |  I_RMS  |  V_MAX  |  V_AVG  |  V_RMS  |  P_swi  |  P_cap  |  P_tra  |  P_tot  |    Pi   |    Po   |   Eta   |')
        print('|----------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|')

        # Primary
        I_pri = [np.max(time['Dc']['i_dc_pri']), np.mean(time['Dc']['i_dc_pri']), rms(time['Dc']['i_dc_pri'])]
        V_pri = [np.max(time['Dc']['v_dc_pri']), np.mean(time['Dc']['v_dc_pri']), rms(time['Dc']['v_dc_pri'])]
        P_pri = [np.sum(P_ALL[0:3, :]), np.sum(P_ALL[8, 0]), np.sum(P_ALL[10, 0]) + np.sum(P_ALL[10, 2])/2]
        Eta_pri = 1 - np.sum(P_pri) / Pin
        print('| PRIMARY  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.0f  |  %5.0f  |  %5.1f  |' % (
            I_pri[0], I_pri[1], I_pri[2], V_pri[0], V_pri[1], V_pri[2], P_pri[0], P_pri[1], P_pri[2], np.sum(P_pri), Pin, Pin-np.sum(P_pri), Eta_pri*100))
        print('|----------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|')

        # Secondary
        I_sec = [np.max(time['Dc']['i_dc_sec']), np.mean(time['Dc']['i_dc_sec']), rms(time['Dc']['i_dc_sec'])]
        V_sec = [np.max(time['Dc']['v_dc_sec']), np.mean(time['Dc']['v_dc_sec']), rms(time['Dc']['v_dc_sec'])]
        P_sec = [np.sum(P_ALL[4:7, :]), np.sum(P_ALL[9, 0]), np.sum(P_ALL[10, 1]) + np.sum(P_ALL[10, 2]) / 2]
        Eta_sec = 1 - np.sum(P_sec)/(I_sec_rms*V_sec_rms)
        print('| SECONDARY|  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.0f  |  %5.0f  |  %5.1f  |' % (
            I_sec[0], I_sec[1], I_sec[2], V_sec[0], V_sec[1], V_sec[2], P_sec[0], P_sec[1], P_sec[2], np.sum(P_sec), Po+np.sum(P_sec), Po, Eta_sec*100))
        print('|----------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|')

        # Total
        Eta_total = 1 - (np.sum(P_pri) + np.sum(P_sec))/Pin
        print(
            '| TOTAL    |         |         |         |         |         |         |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.0f  |  %5.0f  |  %5.1f  |' % (
                P_pri[0]+P_sec[0], P_pri[1]+P_sec[1], P_pri[2]+P_sec[2], np.sum(P_pri) + np.sum(P_sec), Pin, Po, Eta_total * 100))
        print(
            '|----------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|')
        print('\n')

    def _print_table(title, idx_list, cap_idx):
        print(title)
        print(
            '|          |                     CURRENT AND VOLTAGE                   |            AVERAGE LOSSES             |                THERMAL                |   Eta   |')
        print(
            '| item ID  |  I_MAX  |  I_AVG  |  I_RMS  |  V_MAX  |  V_AVG  |  V_RMS  |  P_T_s  |  P_T_c  |  P_D_s  |  P_D_c  | T_T_MAX | T_T_AVG | T_D_MAX | T_D_AVG |    -    |')
        print(
            '|----------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|')

        for j in idx_list:
            print(
                '| %-8s |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.2f  |  %5.2f  |  %5.2f  |  %5.2f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.2f  |' % (
                    lab[j] if j < len(lab) else id1[j], I_ALL[j, 0], I_ALL[j, 1], I_ALL[j, 2], V_ALL[j, 0], V_ALL[j, 1],
                    V_ALL[j, 2], P_ALL[j, 0],
                    P_ALL[j, 1], P_ALL[j, 2], P_ALL[j, 3], T_ALL[j, 0], T_ALL[j, 1], T_ALL[j, 2], T_ALL[j, 3],
                    E_ALL[j, 0] * 100))
            print(
                '|----------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|')

        if cap_idx is not None:
            j = cap_idx
            cap_lab = 'C1' if j == len(id1) else 'C2'
            print(
                '| %-8s |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |                 %5.2f                 |                 %5.1f                 |  %5.2f  |' % (
                    cap_lab, I_ALL[j, 0], I_ALL[j, 1], I_ALL[j, 2], V_ALL[j, 0], V_ALL[j, 1], V_ALL[j, 2],
                    P_ALL[j, 0], T_ALL[j, 0], E_ALL[j, 0] * 100))
            print(
                '|----------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|')
        print('\n')

    _print_overview()
    _print_table("PRIMARY SIDE", list(range(0, 4)), len(id1))
    _print_table("SECONDARY SIDE", list(range(4, 8)), len(id1) + 1)

    ###################################################################################################################
    # MSG Out
    ###################################################################################################################
    print("------------------------------------------")
    print("END: Printing Results")
    print("------------------------------------------")

    ###################################################################################################################
    # Return
    ###################################################################################################################
