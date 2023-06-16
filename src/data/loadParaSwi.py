#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         loadParaSwi
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

# ==============================================================================
# External
# ==============================================================================
import pandas as pd
from os.path import join as pjoin

#######################################################################################################################
# Function
#######################################################################################################################
def loadParaSwi(name, path, setupPara):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("INFO: Loading switch parameters")

    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Variables
    # ==============================================================================
    para = {}
    para['Elec'] = {}
    para['Ther'] = {}
    para['Life'] = {}
    para['Elec']['con'] = {}
    para['Elec']['tab'] = {}
    para['Elec']['vec'] = {}
    para['Ther']['con'] = {}
    para['Ther']['tab'] = {}
    para['Ther']['vec'] = {}

    # ==============================================================================
    # Parameters
    # ==============================================================================#
    col = ['Value-1', 'Value-2', 'Value-3', 'Value-4', 'Value-5', 'Value-6', 'Value-7', 'Value-8', 'Value-9', 'Value-10']
    tab = ['Description', 'Model', 'Symbol', 'Typical', 'Value-1', 'Value-2', 'Value-3', 'Value-4', 'Value-5', 'Value-6']
    lenElec = 24
    lenTher = 6

    ###################################################################################################################
    # Loading Data
    ###################################################################################################################
    # ==============================================================================
    # Path
    # ==============================================================================
    name = name + '.xlsx'
    path = path + '\\Swi'
    filename = pjoin(path, name)

    # ==============================================================================
    # Parameters
    # ==============================================================================
    dataElec = pd.read_excel(filename, sheet_name='electrical')
    dataTher = pd.read_excel(filename, sheet_name='thermal') 

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
    for i in range(0,len(varElecNamesCon)):
        para['Elec']['con'][varElecNamesCon[i]] = varElecValueCon[i]
        para['Elec']['vec'][varElecNamesCon[i]] = varElecValueTab.iloc[i]

    # ------------------------------------------
    # Tabular values
    # ------------------------------------------
    para['Elec']['tab']['Vf'] = dataElec[tab][33:43]
    para['Elec']['tab']['Eon'] = dataElec[tab][53:63]
    para['Elec']['tab']['Eoff'] = dataElec[tab][73:83]
    para['Elec']['tab']['Erec'] = dataElec[tab][93:103]
    para['Elec']['tab']['Vfd'] = dataElec[tab][113:123]
    para['Elec']['tab']['Ciss'] = dataElec[tab][133:143]
    para['Elec']['tab']['Coss'] = dataElec[tab][153:163]
    para['Elec']['tab']['Crss'] = dataElec[tab][173:183]

    # ==============================================================================
    # Thermal
    # ==============================================================================
    for i in range(0,len(varTherNamesCon)):
        para['Ther']['con'][varTherNamesCon[i]] = varTherValueCon[i]
        para['Ther']['vec'][varTherNamesCon[i]] = varTherValueTab.iloc[i]

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
    for i in range(0,len(varElecNamesCon)):
        para['Elec']['vec'][varElecNamesCon[i]] = para['Elec']['vec'][varElecNamesCon[i]].dropna(axis= 0, how='all')

    # Thermal
    for i in range(0,len(varTherNamesCon)):
        para['Ther']['vec'][varTherNamesCon[i]] = para['Ther']['vec'][varTherNamesCon[i]].dropna(axis= 0, how='all')

    # ------------------------------------------
    # Matrix
    # ------------------------------------------
    para['Elec']['tab']['Vf'] = para['Elec']['tab']['Vf'].dropna(axis= 0, how='all')
    para['Elec']['tab']['Vf'] = para['Elec']['tab']['Vf'].dropna(axis= 1, how='all')
    para['Elec']['tab']['Eon'] = para['Elec']['tab']['Eon'].dropna(axis= 0, how='all')
    para['Elec']['tab']['Eon'] = para['Elec']['tab']['Eon'].dropna(axis= 1, how='all')
    para['Elec']['tab']['Eoff'] = para['Elec']['tab']['Eoff'].dropna(axis= 0, how='all')
    para['Elec']['tab']['Eoff'] = para['Elec']['tab']['Eoff'].dropna(axis= 1, how='all')
    para['Elec']['tab']['Erec'] = para['Elec']['tab']['Erec'].dropna(axis= 0, how='all')
    para['Elec']['tab']['Erec'] = para['Elec']['tab']['Erec'].dropna(axis= 1, how='all')
    para['Elec']['tab']['Vfd'] = para['Elec']['tab']['Vfd'].dropna(axis= 0, how='all')
    para['Elec']['tab']['Vfd'] = para['Elec']['tab']['Vfd'].dropna(axis= 1, how='all')
    para['Elec']['tab']['Ciss'] = para['Elec']['tab']['Ciss'].dropna(axis= 0, how='all')
    para['Elec']['tab']['Ciss'] = para['Elec']['tab']['Ciss'].dropna(axis= 1, how='all')
    para['Elec']['tab']['Coss'] = para['Elec']['tab']['Coss'].dropna(axis= 0, how='all')
    para['Elec']['tab']['Coss'] = para['Elec']['tab']['Coss'].dropna(axis= 1, how='all')
    para['Elec']['tab']['Crss'] = para['Elec']['tab']['Crss'].dropna(axis= 0, how='all')
    para['Elec']['tab']['Crss'] = para['Elec']['tab']['Crss'].dropna(axis= 1, how='all')

    # ==============================================================================
    # Losses
    # ==============================================================================
    if setupPara['PWM']['loss'] == 0:
        # ------------------------------------------
        # Matrix
        # ------------------------------------------
        para['Elec']['tab']['Vf'] = para['Elec']['tab']['Vf']*0
        para['Elec']['tab']['Eon'] = para['Elec']['tab']['Eon']*0
        para['Elec']['tab']['Eoff'] = para['Elec']['tab']['Eoff']*0
        para['Elec']['tab']['Erec'] = para['Elec']['tab']['Erec']*0
        para['Elec']['tab']['Vfd'] = para['Elec']['tab']['Vfd']*0
        para['Elec']['tab']['Ciss'] = para['Elec']['tab']['Ciss']*0
        para['Elec']['tab']['Coss'] = para['Elec']['tab']['Coss']*0
        para['Elec']['tab']['Crss'] = para['Elec']['tab']['Crss']*0

        # ------------------------------------------
        # Constant
        # ------------------------------------------
        para['Elec']['con']['Vf'] = para['Elec']['con']['Vf']*0
        para['Elec']['con']['Eon'] = para['Elec']['con']['Eon']*0
        para['Elec']['con']['Eoff'] = para['Elec']['con']['Eoff']*0
        para['Elec']['con']['Erec'] = para['Elec']['con']['Erec']*0
        para['Elec']['con']['Vfd'] = para['Elec']['con']['Vfd']*0
        para['Elec']['con']['Ron'] = para['Elec']['con']['Ron']*0
        para['Elec']['con']['Roff'] = para['Elec']['con']['Roff']*0
        para['Elec']['con']['RonD'] = para['Elec']['con']['RonD']*0
        para['Elec']['con']['RoffD'] = para['Elec']['con']['RoffD']*0
        para['Elec']['con']['Ciss'] = para['Elec']['con']['Ciss']*0
        para['Elec']['con']['Coss'] = para['Elec']['con']['Coss']*0
        para['Elec']['con']['Crss'] = para['Elec']['con']['Crss']*0

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return para