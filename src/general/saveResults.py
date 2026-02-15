#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         saveResults
# Date:         01.05.2024
# Author:       Dr. Pascal A. Schirmer
# Version:      V.1.0
# Copyright:    Pascal Schirmer
#######################################################################################################################
#######################################################################################################################

#######################################################################################################################
# Function Description
#######################################################################################################################
"""
This function saves the results for the different simulation setups and topologies.
Inputs:     1) time:    time domain results
            2) freq:    frequency domain results
            3) sweep:   sweep domain results 
            4) setup:   includes all setup variables 
            5) path:    path variables
"""

#######################################################################################################################
# Import libs
#######################################################################################################################
# ==============================================================================
# Internal
# ==============================================================================

# ==============================================================================
# External
# ==============================================================================
import os
import pandas as pd
import numpy as np
from datetime import datetime
from pathlib import Path
import json


#######################################################################################################################
# Function
#######################################################################################################################
def saveResults(time, freq, sweep, setup, path):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("------------------------------------------")
    print("START: Saving Results")
    print("------------------------------------------")

    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Labels
    # ==============================================================================
    if setup['Top']['sourceType'] == 'B2':
        idSw = ['S1', 'S2']
    elif setup['Top']['sourceType'] == 'B4':
        idSw = ['S1', 'S2', 'S3', 'S4']
    elif setup['Top']['sourceType'] == 'B6':
        idSw = ['S1', 'S2', 'S3', 'S4', 'S5', 'S6']
    elif setup['Top']['sourceType'] == 'DAB':
        idSw = ['S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8']
    else:
        print("WARN: Invalid topology assuming B6")
        idSw = ['S1', 'S2', 'S3', 'S4', 'S5', 'S6']

    # ==============================================================================
    # General
    # ==============================================================================
    now = datetime.now()
    dir_name = setup['Exp']['Name']
    dt_string = now.strftime("%d%m%Y%H%M%S")

    # ==============================================================================
    # Variables
    # ==============================================================================
    fel = setup['Top']['fel']
    fsim = setup['Exp']['fsim']
    f = fsim * np.linspace(0, 0.5, int(len(time['Sw']['t']) / 2)) / fel
    setup['Dat']['stat']['Io'] = abs(setup['Dat']['stat']['Io'])
    setup['Dat']['stat']['So'] = abs(setup['Dat']['stat']['So'])

    ###################################################################################################################
    # Directory
    ###################################################################################################################
    path = path['resPath'] + '\\' + dir_name
    Path(path).mkdir(parents=True, exist_ok=True)
    print("INFO: Generating directory")

    ####################################################################################################################
    # Save Setup
    ####################################################################################################################
    setup_name = 'setup_' + dir_name + '_' + dt_string + '.txt'
    os.chdir(path)
    with open(setup_name, 'w') as file:
        file.write(json.dumps(setup['Exp']))
        file.write(json.dumps(setup['Dat']))
        file.write(json.dumps(setup['Par']))
        file.write(json.dumps(setup['Top']))
    print("INFO: Saving setup files")

    ####################################################################################################################
    # Output
    ####################################################################################################################
    # ==============================================================================
    # Time Data
    # ==============================================================================
    try:
        # ------------------------------------------
        # Phases
        # ------------------------------------------
        try:
            nameTime = 'time_' + dir_name + '_' + dt_string + '_AcDc.xlsx'
            with pd.ExcelWriter(nameTime) as writer:
                time['Sw']['t'].to_excel(writer, sheet_name='time')
                time['Sw'].to_excel(writer, sheet_name='Sw')
                pd.DataFrame.from_dict(time['Ac']).to_excel(writer, sheet_name='Ac')
                pd.DataFrame.from_dict(time['Dc']).to_excel(writer, sheet_name='Dc')
            print("INFO: Saving per phase results")
        except:
            print("WARN: Saving per phase results failed")

        # ------------------------------------------
        # Switches
        # ------------------------------------------
        try:
            for i in range(0, len(idSw)):
                nameTime = 'time_' + dir_name + '_' + dt_string + '_' + idSw[i] + '.xlsx'
                with pd.ExcelWriter(nameTime) as writer:
                    time['Sw']['t'].to_excel(writer, sheet_name='time')
                    time['Elec']['sw'][idSw[i]].to_excel(writer, sheet_name='elec')
                    time['Loss']['sw'][idSw[i]].to_excel(writer, sheet_name='loss')
                    pd.DataFrame.from_dict(time['Ther']['sw']).to_excel(writer, sheet_name='ther')
            print("INFO: Saving per switch results")
        except:
            print("WARN: Saving per switch results failed")

        # ------------------------------------------
        # Capacitor
        # ------------------------------------------
        try:
            nameTime = 'time_' + dir_name + '_' + dt_string + '_C1.xlsx'
            with pd.ExcelWriter(nameTime) as writer:
                time['Sw']['t'].to_excel(writer, sheet_name='time')
                time['Elec']['cap']['C1'].to_excel(writer, sheet_name='elec')
                time['Loss']['cap']['C1'].to_excel(writer, sheet_name='loss')
                pd.DataFrame.from_dict(time['Ther']['cap']['C1']).to_excel(writer, sheet_name='ther')
            print("INFO: Saving per capacitor results")
        except:
            print("WARN: Saving per capacitor results failed")

        # ------------------------------------------
        # Output
        # ------------------------------------------
        print("INFO: Time domain data saved")

    except:
        print("WARN: Couldn't save time domain data")

    # ==============================================================================
    # Freq Data
    # ==============================================================================
    try:
        # ------------------------------------------
        # Output
        # ------------------------------------------
        nameFreq = 'freq_' + dir_name + '_' + dt_string + '.xlsx'
        with pd.ExcelWriter(nameFreq) as writer:
            pd.DataFrame.from_dict(f).to_excel(writer, sheet_name='freq')
            pd.DataFrame.from_dict(freq['Sw']).to_excel(writer, sheet_name='Sw')
            pd.DataFrame.from_dict(freq['Ac']).to_excel(writer, sheet_name='Ac')
            pd.DataFrame.from_dict(freq['Dc']).to_excel(writer, sheet_name='Dc')

        # ------------------------------------------
        # Output
        # ------------------------------------------
        print("INFO: Frequency domain data saved")
    except:
        print("WARN: Couldn't save frequency domain data")

    # ==============================================================================
    # Sweep Data
    # ==============================================================================
    try:
        # ------------------------------------------
        # Output
        # ------------------------------------------
        nameSweep = 'sweep_' + dir_name + '_' + dt_string + '.xlsx'
        with pd.ExcelWriter(nameSweep) as writer:
            pd.DataFrame.from_dict(sweep['Mi']).to_excel(writer, sheet_name='Mi')
            pd.DataFrame.from_dict(sweep['Ac']['ana']).to_excel(writer, sheet_name='Ac_ana')
            pd.DataFrame.from_dict(sweep['Ac']['num']).to_excel(writer, sheet_name='Ac_num')
            pd.DataFrame.from_dict(sweep['Dc']['ana']).to_excel(writer, sheet_name='Dc_ana')
            pd.DataFrame.from_dict(sweep['Dc']['num']).to_excel(writer, sheet_name='Dc_num')

        # ------------------------------------------
        # Output
        # ------------------------------------------
        print("INFO: Sweeping domain data saved")
    except:
        print("WARN: Couldn't save sweeping domain data")

    ###################################################################################################################
    # MSG Out
    ###################################################################################################################
    print("------------------------------------------")
    print("END: Saving Results")
    print("------------------------------------------")

    ###################################################################################################################
    # Return
    ###################################################################################################################
