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
import math


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
        if setup['Top']['sourceType'] == 'B2':
            id1 = ['S1', 'S2']
        elif setup['Top']['sourceType'] == 'B4':
            id1 = ['S1', 'S2', 'S3', 'S4']
        elif setup['Top']['sourceType'] == 'B6':
            id1 = ['S1', 'S2', 'S3', 'S4', 'S5', 'S6']
        elif setup['Top']['sourceType'] == 'DAB':
            id1 = ['S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8']
        else:
            print("WARN: Invalid topology assuming B6")
            id1 = ['S1', 'S2', 'S3', 'S4', 'S5', 'S6']

    id2 = ['T' + s[1:] for s in id1]
    id3 = ['D' + s[1:] for s in id1]
    lab = id1 + ['C1']

    # ==============================================================================
    # Variables
    # ==============================================================================
    pt_sec = time['Dc']['i_dc_sec'] * time['Dc']['v_dc_sec'] 
    pt_pri = time['Dc']['i_dc_pri'] * time['Dc']['v_dc_pri']

    # ==============================================================================
    # Init
    # ==============================================================================
    I_ALL = np.zeros((len(id1) + 1, 3))
    V_ALL = np.zeros((len(id1) + 1, 3))
    P_ALL = np.zeros((len(id1) + 1, 4))
    T_ALL = np.zeros((len(id1) + 1, 4))
    E_ALL = np.zeros((len(id1) + 1, 1))

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

        # Currents
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
    I_ALL[-1, 0] = np.max(time['Elec']['cap']['C1']['i_c'])
    I_ALL[-1, 1] = np.mean(time['Elec']['cap']['C1']['i_c'])
    I_ALL[-1, 2] = rms(time['Elec']['cap']['C1']['i_c'])

    # Currents
    V_ALL[-1, 0] = np.max(time['Elec']['cap']['C1']['v_c'])
    V_ALL[-1, 1] = np.mean(time['Elec']['cap']['C1']['v_c'])
    V_ALL[-1, 2] = rms(time['Elec']['cap']['C1']['v_c'])

    # Losses
    P_ALL[-1, 0] = np.mean(time['Loss']['cap']['C1']['p_L']) / 4
    P_ALL[-1, 1] = np.mean(time['Loss']['cap']['C1']['p_L']) / 4
    P_ALL[-1, 2] = np.mean(time['Loss']['cap']['C1']['p_L']) / 4
    P_ALL[-1, 3] = np.mean(time['Loss']['cap']['C1']['p_L']) / 4

    # Thermal
    T_ALL[-1, 0] = np.max(time['Ther']['cap']['C1'])
    T_ALL[-1, 1] = np.mean(time['Ther']['cap']['C1'])
    T_ALL[-1, 2] = np.max(time['Ther']['cap']['C1'])
    T_ALL[-1, 3] = np.mean(time['Ther']['cap']['C1'])

    # Efficiency
    E_ALL[-1, 0] = np.abs(rms(pt_sec) - np.mean(time['Loss']['cap']['C1']['p_L'])) / rms(pt_sec)

    # ==============================================================================
    # Converter Level
    # ==============================================================================
    # ------------------------------------------
    # Currents and Voltages
    # ------------------------------------------
    I_pri_rms = rms(time['Ac']['i_ac_pri'])
    I_sec_rms = rms(time['Ac']['i_ac_sec'])
    V_pri_rms = rms(time['Ac']['v_ac_pri'])
    V_sec_rms = rms(time['Ac']['v_ac_sec'])
    I_dc_in_rms = rms(time['Dc']['i_dc_pri'])
    I_dc_out_rms = rms(time['Dc']['i_dc_sec'])
    V_dc_in_rms = rms(time['Dc']['v_dc_pri'])
    V_dc_out_rms = rms(time['Dc']['v_dc_sec'])

    # ------------------------------------------
    # Power
    # ------------------------------------------
    Pin = I_dc_in_rms * V_dc_in_rms
    Po = I_dc_out_rms * V_dc_out_rms
    Eta = Po / Pin if Pin != 0 else 0

    # ------------------------------------------
    # Others
    # ------------------------------------------
    # Currents
    I_CON = [np.max(I_ALL[:, 0]), np.mean(I_ALL[:, 1]), rms(I_ALL[:, 2])]

    # Voltages
    V_CON = [np.max(V_ALL[:, 0]), np.mean(V_ALL[:, 1]), rms(V_ALL[:, 2])]

    # Losses
    P_CON = [np.sum(P_ALL[:, 0]), np.sum(P_ALL[:, 1]), np.sum(P_ALL[:, 2]), np.sum(P_ALL[:, 3])]

    # Thermal
    T_CON = [np.max(T_ALL[:, 0]), np.mean(T_ALL[:, 1]), np.max(T_ALL[:, 2]), np.mean(T_ALL[:, 3])]

    # Efficiency
    E_CON = np.abs(Pin - np.sum(np.sum(P_ALL))) / Pin

    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    def _print_table(title, idx_list, include_cap, pt_ref, meta=None):
        if meta is not None:
            print(title + " META (RMS)")
            print(f"| {meta['label']:10s} | I_dc={meta['I_dc']:.3f} A | V_dc={meta['V_dc']:.3f} V | I_ac={meta['I_ac']:.3f} A | V_ac={meta['V_ac']:.3f} V | P={meta['P']:.3f} W | Eta={meta['Eta']:.3f} % |")
            print('|--------------------------------------------------------------------------|')
        print(title)
        print(
            '|          |                     CURRENT AND VOLTAGE                   |            AVERAGE LOSSES             |                THERMAL                |   Eta   |')
        print(
            '| item ID  |  I_MAX  |  I_AVG  |  I_RMS  |  V_MAX  |  V_AVG  |  V_RMS  |  P_T_s  |  P_T_c  |  P_D_s  |  P_D_c  | T_T_MAX | T_T_AVG | T_D_MAX | T_D_AVG |    -    |')
        print(
            '|----------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|')

        for i in idx_list:
            print(
                '| %-8s |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.2f  |  %5.2f  |  %5.2f  |  %5.2f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.2f  |' % (
                    lab[i], I_ALL[i, 0], I_ALL[i, 1], I_ALL[i, 2], V_ALL[i, 0], V_ALL[i, 1], V_ALL[i, 2], P_ALL[i, 0],
                    P_ALL[i, 1], P_ALL[i, 2], P_ALL[i, 3], T_ALL[i, 0], T_ALL[i, 1], T_ALL[i, 2], T_ALL[i, 3],
                    E_ALL[i, 0] * 100))
            print(
                '|----------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|')

        if include_cap:
            i = len(id1)
            print(
                '| %-8s |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |                 %5.2f                 |                 %5.1f                 |  %5.2f  |' % (
                    lab[i], I_ALL[i, 0], I_ALL[i, 1], I_ALL[i, 2], V_ALL[i, 0], V_ALL[i, 1], V_ALL[i, 2],
                    P_ALL[i, 0] * 4, T_ALL[i, 0], E_ALL[i, 0] * 100))
            print(
                '|----------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|')

        # Totals
        I_CON = [np.max(I_ALL[idx_list, 0]), np.mean(I_ALL[idx_list, 1]), rms(I_ALL[idx_list, 2])]
        V_CON = [np.max(V_ALL[idx_list, 0]), np.mean(V_ALL[idx_list, 1]), rms(V_ALL[idx_list, 2])]
        P_CON = [np.sum(P_ALL[idx_list, 0]), np.sum(P_ALL[idx_list, 1]), np.sum(P_ALL[idx_list, 2]), np.sum(P_ALL[idx_list, 3])]
        T_CON = [np.max(T_ALL[idx_list, 0]), np.mean(T_ALL[idx_list, 1]), np.max(T_ALL[idx_list, 2]), np.mean(T_ALL[idx_list, 3])]
        if pt_ref is not None and np.mean(pt_ref) != 0:
            P_out = np.abs(np.mean(pt_ref))
            P_loss = np.sum(np.sum(P_ALL[idx_list, :]))
            if include_cap:
                P_loss = P_loss + np.mean(time['Loss']['cap']['C1']['p_L'])
            E_CON = (P_out - P_loss) / P_out
        else:
            E_CON = 0

        print(
            '|  Total   |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.1f  |  %5.2f  |' % (
                I_CON[0], I_CON[1], I_CON[2], V_CON[0], V_CON[1], V_CON[2], P_CON[0], P_CON[1], P_CON[2], P_CON[3],
                T_CON[0], T_CON[1], T_CON[2], T_CON[3], E_CON * 100))
        print(
            '|----------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|')


    meta_pri = {
        'label': 'PRIMARY',
        'I_dc': I_dc_in_rms,
        'V_dc': V_dc_in_rms,
        'I_ac': I_pri_rms,
        'V_ac': V_pri_rms,
        'P': Pin,
        'Eta': 100 * Po / Pin if Pin != 0 else 0
    }
    meta_sec = {
        'label': 'SECONDARY',
        'I_dc': I_dc_out_rms,
        'V_dc': V_dc_out_rms,
        'I_ac': I_sec_rms,
        'V_ac': V_sec_rms,
        'P': Po,
        'Eta': 100 * Po / Pin if Pin != 0 else 0
    }

    _print_table("PRIMARY SIDE", list(range(0, 4)), False, pt_pri, meta=meta_pri)
    _print_table("SECONDARY SIDE", list(range(4, 8)), True, pt_sec, meta=meta_sec)

    ###################################################################################################################
    # MSG Out
    ###################################################################################################################
    print("------------------------------------------")
    print("END: Printing Results")
    print("------------------------------------------")

    ###################################################################################################################
    # Return
    ###################################################################################################################
