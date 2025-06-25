#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         loadParaTra
# Date:         08.05.2025
# Author:       Max V. Mueller
# Version:      V.1.0
# Copyright:    Pascal Schirmer
#######################################################################################################################
#######################################################################################################################

#######################################################################################################################
# Function Description
#######################################################################################################################
"""
This function loads the parameters for the transformer under \par. This includes experimental, data, topology, and
electrical as well as thermal parameter information. The parameters are summarized in one common para variable.
Inputs:     1) name:    name of the parameter file for the transformer
            2) path:    includes all path variables
            2) setup:   includes all simulation variables
Outputs:    1) para:    output parameter file for the transformer
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
import numpy as np
from os.path import join as pjoin


#######################################################################################################################
# Function
#######################################################################################################################
def loadParaTra(name, path, setup):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("INFO: Loading transformer parameters")

    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Variables
    # ==============================================================================
    para = {'Elec': {}, 'Ther': {}, 'Life': {}, 'SS':  {}}
    para['Elec']['con'] = {}
    para['Elec']['vec'] = {}
    para['Elec']['tab'] = {}
    para['Ther']['con'] = {}
    para['Ther']['vec'] = {}
    para['Ther']['tab'] = {}
    para['Life']['con'] = {}
    para['Life']['tab'] = {}

    # ==============================================================================
    # Parameters
    # ==============================================================================
    col = ['Value-1', 'Value-2', 'Value-3', 'Value-4', 'Value-5', 'Value-6', 'Value-7', 'Value-8', 'Value-9',
           'Value-10']
    tab = ['Description', 'Model', 'Symbol', 'Typical', 'Value-1', 'Value-2', 'Value-3', 'Value-4', 'Value-5',
           'Value-6']
    lenElec = 22
    lenTher = 12

    ###################################################################################################################
    # Loading Data
    ###################################################################################################################
    # ==============================================================================
    # Path
    # ==============================================================================
    name = name + '.xlsx'
    filename = pjoin(path, 'Tra', name)

    # ==============================================================================
    # Parameters
    # ==============================================================================
    dataElec = pd.read_excel(filename, sheet_name='electrical')
    dataTher = pd.read_excel(filename, sheet_name='thermal')
    dataSS = pd.read_excel(filename, sheet_name='statespace')

    ###################################################################################################################
    # Pre-Processing
    ###################################################################################################################
    # ==============================================================================
    # Electrical
    # ==============================================================================
    varElecNamesCon = dataElec['Symbol'][0:lenElec]
    varElecValueCon = dataElec['Typical'][0:lenElec]
    varElecValueTab = dataElec.loc[:, col][0:lenElec]

    # ==============================================================================
    # Thermal
    # ==============================================================================
    varTherNamesCon = dataTher['Symbol'][0:lenTher]
    varTherValueCon = dataTher['Typical'][0:lenTher]
    varTherValueTab = dataTher.loc[:, col][0:lenTher]

    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    # ==============================================================================
    # Electrical
    # ==============================================================================
    # ------------------------------------------
    # Constant values
    # ------------------------------------------
    for i in range(0, len(varElecNamesCon)):
        para['Elec']['con'][varElecNamesCon[i]] = varElecValueCon[i]
        para['Elec']['vec'][varElecNamesCon[i]] = varElecValueTab.iloc[i]

    # ------------------------------------------
    # Tabular values
    # ------------------------------------------
    para['Elec']['tab']['alpha'] = dataElec[tab][31:41]
    para['Elec']['tab']['beta'] = dataElec[tab][51:61]
    para['Elec']['tab']['Pvsin'] = dataElec[tab][71:81]
    #para['Elec']['tab']['tbd'] = dataElec[tab][91:101]

    # ==============================================================================
    # Thermal
    # ==============================================================================
    for i in range(0, len(varTherNamesCon)):
        para['Ther']['con'][varTherNamesCon[i]] = varTherValueCon[i]
        para['Ther']['vec'][varTherNamesCon[i]] = varTherValueTab.iloc[i]

    # ==============================================================================
    # State Space
    # ==============================================================================
    # Ports:
    para['SS']['PORT'] = {}

    # Matrix A
    A = dataSS.iloc[0:8, 1:9]
    A = A.dropna(how='all').reset_index(drop=True)
    A = A.dropna(axis=1, how='all')
    para['SS']['PORT']['A'] = A.to_numpy()

    # Vector B
    B = dataSS.iloc[9:17, 1]
    B = B.dropna(how='all').reset_index(drop=True)
    B = B.to_numpy()
    para['SS']['PORT']['B'] = B[:, np.newaxis]

    # Matrix C
    C = dataSS.iloc[18:21, 1:9]
    C = C.dropna(how='all').reset_index(drop=True)
    C = C.dropna(axis=1, how='all')
    para['SS']['PORT']['C'] = C.to_numpy()

    # Vector D
    D = dataSS.iloc[22:25, 1]
    D = D.dropna(how='all').reset_index(drop=True)
    D = D.to_numpy()
    para['SS']['PORT']['D'] = D[:, np.newaxis]

    # Winding values:
    para['SS']['WDG'] = {}

    # Matrix A
    A = dataSS.iloc[0:8, 15:23]
    A = A.dropna(how='all').reset_index(drop=True)
    A = A.dropna(axis=1, how='all')
    para['SS']['WDG']['A'] = A.to_numpy()

    # Vector B
    B = dataSS.iloc[9:17, 15]
    B = B.dropna(how='all').reset_index(drop=True)
    B = B.to_numpy()
    para['SS']['WDG']['B'] = B[:, np.newaxis]

    # Matrix C
    C = dataSS.iloc[18:21, 15:23]
    C = C.dropna(how='all').reset_index(drop=True)
    C = C.dropna(axis=1, how='all')
    para['SS']['WDG']['C'] = C.to_numpy()

    # Vector D
    D = dataSS.iloc[22:25, 15]
    D = D.dropna(how='all').reset_index(drop=True)
    D = D.to_numpy()
    para['SS']['WDG']['D'] = D[:, np.newaxis]
    ###################################################################################################################
    # Post-Processing
    ###################################################################################################################
    # ==============================================================================
    # Remove NaNs
    # ==============================================================================
    # ------------------------------------------
    # Vectors
    # ------------------------------------------
    # Electrical
    for i in range(0, len(varElecNamesCon)):
        para['Elec']['vec'][varElecNamesCon[i]] = para['Elec']['vec'][varElecNamesCon[i]].dropna(axis=0, how='all')

    # Thermal
    for i in range(0, len(varTherNamesCon)):
        para['Ther']['vec'][varTherNamesCon[i]] = para['Ther']['vec'][varTherNamesCon[i]].dropna(axis=0, how='all')

    # ------------------------------------------
    # Matrix
    # ------------------------------------------
    para['Elec']['tab']['alpha'] = para['Elec']['tab']['alpha'].dropna(axis=0, how='all')
    para['Elec']['tab']['alpha'] = para['Elec']['tab']['alpha'].dropna(axis=1, how='all')
    para['Elec']['tab']['beta'] = para['Elec']['tab']['beta'].dropna(axis=0, how='all')
    para['Elec']['tab']['beta'] = para['Elec']['tab']['beta'].dropna(axis=1, how='all')
    para['Elec']['tab']['Pvsin'] = para['Elec']['tab']['Pvsin'].dropna(axis=0, how='all')
    para['Elec']['tab']['Pvsin'] = para['Elec']['tab']['Pvsin'].dropna(axis=1, how='all')

    """ 
    ToDo
    # ==============================================================================
    # Losses
    # ==============================================================================
    if setup['Par']['PWM']['loss'] == 0:
        # ------------------------------------------
        # Matrix
        # ------------------------------------------
        para['Elec']['tab']['ESR'] = para['Elec']['tab']['ESR'] * 0
        para['Elec']['tab']['tan'] = para['Elec']['tab']['tan'] * 0
        para['Elec']['tab']['Kr'] = para['Elec']['tab']['Kr'] * 0

        # ------------------------------------------
        # Constant
        # ------------------------------------------
        para['Elec']['con']['ESR'] = para['Elec']['con']['ESR'] * 0
        para['Elec']['con']['tan'] = para['Elec']['con']['tan'] * 0
        para['Elec']['con']['Kr'] = para['Elec']['con']['Kr'] * 0
    """

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return para
