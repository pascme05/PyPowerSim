#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         loadSetup
# Date:         27.04.2024
# Author:       Dr. Pascal A. Schirmer
# Version:      V.1.0
# Copyright:    Pascal Schirmer
#######################################################################################################################
#######################################################################################################################

#######################################################################################################################
# Function Description
#######################################################################################################################
"""
This function loads the remaining setup parameters from the setup file located under \config. This includes experimental,
data, topology, and electrical as well as thermal parameter information. The parameters are summarized in one common
setup variable.
Inputs:     1) setup:   includes all simulation variables
            2) path:    includes all path variables
Outputs:    1) setup:   extended setup variable
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
import pandas as pd
from os.path import join as pjoin


#######################################################################################################################
# Function
#######################################################################################################################
def loadSetup(setup, path):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("------------------------------------------")
    print("Loading Setup")
    print("------------------------------------------")
    print("START: Loading setup")

    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    column = 'Value'

    ###################################################################################################################
    # Loading
    ###################################################################################################################
    # ==============================================================================
    # Path and Filename
    # ==============================================================================
    name = setup['Exp']['conf'] + '.xlsx'
    path = path['conPath']
    filename = pjoin(path, name)

    # ==============================================================================
    # Loading Config
    # ==============================================================================
    try:
        setupExpRaw = pd.read_excel(filename, sheet_name='Exp')
        setupDatRaw = pd.read_excel(filename, sheet_name='Dat')
        setupTopRaw = pd.read_excel(filename, sheet_name='Top')
        setupParRaw = pd.read_excel(filename, sheet_name='Par')
        setupMagRaw = pd.read_excel(filename, sheet_name='Mag')
        print("INFO: Setup file loaded")
    except:
        setupExpRaw = []
        setupDatRaw = []
        setupTopRaw = []
        setupParRaw = []
        setupMagRaw = []
        print("ERROR: Setup file could not be loaded")

    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    # ==============================================================================
    # Experiment struct
    # ==============================================================================
    try:
        for i in range(0, setupExpRaw.shape[0]):
            setup['Exp'][setupExpRaw['Variable'][i]] = setupExpRaw[column][i]
        print("INFO: Experimental setup file loaded")
    except:
        print("ERROR: Experimental setup file could not be loaded")

    # ==============================================================================
    # Data struct
    # ==============================================================================
    try:
        for i in range(0, setupDatRaw.shape[0]):
            if setupDatRaw['Category'][i] == 'stat':
                setup['Dat']['stat'][setupDatRaw['Variable'][i]] = setupDatRaw[column][i]
            else:
                setup['Dat']['trans'][setupDatRaw['Variable'][i]] = setupDatRaw[column][i]
        print("INFO: Data setup file loaded")
    except:
        print("ERROR: Data setup file could not be loaded")

    # ==============================================================================
    # Topology struct
    # ==============================================================================
    try:
        for i in range(0, setupTopRaw.shape[0]):
            setup['Top'][setupTopRaw['Variable'][i]] = setupTopRaw[column][i]
        print("INFO: Topology setup file loaded")
    except:
        print("ERROR: Topology setup file could not be loaded")

    # ==============================================================================
    # Parameter struct
    # ==============================================================================
    try:
        for i in range(0, setupParRaw.shape[0]):
            setup['Par'][setupParRaw['Category'][i]][setupParRaw['Variable'][i]] = setupParRaw[column][i]
        print("INFO: Parameter setup file loaded")
    except:
        print("ERROR: Parameter setup file could not be loaded")

    # ==============================================================================
    # Magnetics topology (short circuit/open circuit/load): Add to Topology struct
    # ==============================================================================
    try:
        for i in range(0, setupMagRaw.shape[0]):
            setup['Top'][setupMagRaw['Variable'][i]] = setupMagRaw[column][i]
        print("INFO: Parameter setup file loaded")
    except:
        print("ERROR: Parameter setup file could not be loaded")

    ###################################################################################################################
    # MSG Out
    ###################################################################################################################
    print("DONE: Loading setup")
    print("\n")

    return setup
