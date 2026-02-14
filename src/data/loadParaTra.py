#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         loadParaTra
# Date:         07.02.2026
# Author:       Dr. Pascal A. Schirmer
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
from scipy.interpolate import griddata, RegularGridInterpolator


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
    para = {'Elec': {}, 'Ther': {}, 'Life': {}}
    para['Elec']['con'] = {}
    para['Elec']['tab'] = {}
    para['Elec']['vec'] = {}
    para['Ther']['con'] = {}
    para['Ther']['tab'] = {}
    para['Ther']['vec'] = {}

    # ==============================================================================
    # Parameters
    # ==============================================================================#
    col = ['Value-1', 'Value-2', 'Value-3', 'Value-4', 'Value-5', 'Value-6', 'Value-7', 'Value-8', 'Value-9',
           'Value-10']
    tab = ['Description', 'Model', 'Symbol', 'Typical', 'Value-1', 'Value-2', 'Value-3', 'Value-4', 'Value-5',
           'Value-6']
    lenElec = 23
    lenTher = 8

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
    para['Elec']['tab']['Rp'] = dataElec[tab][33:43]
    para['Elec']['tab']['Rs'] = dataElec[tab][53:63]
    para['Elec']['tab']['Lm'] = dataElec[tab][73:83]
    para['Elec']['tab']['Rc1'] = dataElec[tab][93:103]
    para['Elec']['tab']['Rc2'] = dataElec[tab][113:123]
    para['Elec']['tab']['Cp'] = dataElec[tab][133:143]
    para['Elec']['tab']['Cs'] = dataElec[tab][153:163]
    para['Elec']['tab']['Cps'] = dataElec[tab][173:183]

    # ==============================================================================
    # Thermal
    # ==============================================================================
    for i in range(0, len(varTherNamesCon)):
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
    for i in range(0, len(varElecNamesCon)):
        para['Elec']['vec'][varElecNamesCon[i]] = para['Elec']['vec'][varElecNamesCon[i]].dropna(axis=0, how='all')

    # Thermal
    for i in range(0, len(varTherNamesCon)):
        para['Ther']['vec'][varTherNamesCon[i]] = para['Ther']['vec'][varTherNamesCon[i]].dropna(axis=0, how='all')

    # ------------------------------------------
    # Matrix
    # ------------------------------------------
    para['Elec']['tab']['Rp'] = para['Elec']['tab']['Rp'].dropna(axis=0, how='all')
    para['Elec']['tab']['Rp'] = para['Elec']['tab']['Rp'].dropna(axis=1, how='all')
    para['Elec']['tab']['Rs'] = para['Elec']['tab']['Rs'].dropna(axis=0, how='all')
    para['Elec']['tab']['Rs'] = para['Elec']['tab']['Rs'].dropna(axis=1, how='all')
    para['Elec']['tab']['Lm'] = para['Elec']['tab']['Lm'].dropna(axis=0, how='all')
    para['Elec']['tab']['Lm'] = para['Elec']['tab']['Lm'].dropna(axis=1, how='all')
    para['Elec']['tab']['Rc1'] = para['Elec']['tab']['Rc1'].dropna(axis=0, how='all')
    para['Elec']['tab']['Rc1'] = para['Elec']['tab']['Rc1'].dropna(axis=1, how='all')
    para['Elec']['tab']['Rc2'] = para['Elec']['tab']['Rc2'].dropna(axis=0, how='all')
    para['Elec']['tab']['Rc2'] = para['Elec']['tab']['Rc2'].dropna(axis=1, how='all')
    para['Elec']['tab']['Cp'] = para['Elec']['tab']['Cp'].dropna(axis=0, how='all')
    para['Elec']['tab']['Cp'] = para['Elec']['tab']['Cp'].dropna(axis=1, how='all')
    para['Elec']['tab']['Cs'] = para['Elec']['tab']['Cs'].dropna(axis=0, how='all')
    para['Elec']['tab']['Cs'] = para['Elec']['tab']['Cs'].dropna(axis=1, how='all')
    para['Elec']['tab']['Cps'] = para['Elec']['tab']['Cps'].dropna(axis=0, how='all')
    para['Elec']['tab']['Cps'] = para['Elec']['tab']['Cps'].dropna(axis=1, how='all')

    # ==============================================================================
    # Losses
    # ==============================================================================
    if setup['Par']['PWM']['loss'] == 0:
        # ------------------------------------------
        # Matrix
        # ------------------------------------------
        para['Elec']['tab']['Rp']  = para['Elec']['tab']['Rp'] * 0
        para['Elec']['tab']['Rs']  = para['Elec']['tab']['Rs'] * 0
        para['Elec']['tab']['Rc1'] = para['Elec']['tab']['Rc1'] * 0
        para['Elec']['tab']['Rc2'] = para['Elec']['tab']['Rc2'] * 0
        para['Elec']['tab']['Cp']  = para['Elec']['tab']['Cp'] * 0
        para['Elec']['tab']['Cs']  = para['Elec']['tab']['Cs'] * 0
        para['Elec']['tab']['Cps'] = para['Elec']['tab']['Cps'] * 0

        # ------------------------------------------
        # Constant
        # ------------------------------------------
        para['Elec']['con']['Rp']  = para['Elec']['con']['Rp'] * 0
        para['Elec']['con']['Rs']  = para['Elec']['con']['Rs'] * 0
        para['Elec']['con']['Rc1'] = para['Elec']['con']['Rc1'] * 0
        para['Elec']['con']['Rc2'] = para['Elec']['con']['Rc2'] * 0
        para['Elec']['con']['Cp']  = para['Elec']['con']['Cp'] * 0
        para['Elec']['con']['Cs']  = para['Elec']['con']['Cs'] * 0
        para['Elec']['con']['Cps'] = para['Elec']['con']['Cps'] * 0

    # ==============================================================================
    # Interpolation Functions
    # ==============================================================================
    try:
        # ------------------------------------------
        # Init
        # ------------------------------------------
        nInt = setup['Exp']['int']

        # ------------------------------------------
        # Electrical (Frequency dependency)
        # ------------------------------------------
        # Init
        x = (para['Elec']['vec']['Tj'].to_numpy() * np.ones(
            (len(para['Elec']['vec']['f'].to_numpy()), len(para['Elec']['vec']['Tj'].to_numpy())))).flatten(order='F')
        y = np.tile(para['Elec']['vec']['f'].to_numpy(), len(para['Elec']['vec']['Tj'].to_numpy()))
        xi = np.linspace(np.min(x), np.max(x), nInt)
        yi = np.linspace(np.min(y), np.max(y), nInt)
        xg, yg = np.meshgrid(yi, xi)

        # Primary Winding
        try:
            z = para['Elec']['tab']['Rp'].to_numpy('float').flatten(order='F')
            zi = griddata((y, x), z, (xg, yg), method='linear')
        except:
            zi = np.zeros((len(xi), len(yi)))
        para['Elec']['tab']['Rp_2d'] = RegularGridInterpolator((xi, yi), zi, bounds_error=False, fill_value=None)

        # Secondary Winding
        try:
            z = para['Elec']['tab']['Rs'].to_numpy('float').flatten(order='F')
            zi = griddata((y, x), z, (xg, yg), method='linear')
        except:
            zi = np.zeros((len(xi), len(yi)))
        para['Elec']['tab']['Rs_2d'] = RegularGridInterpolator((xi, yi), zi, bounds_error=False, fill_value=None)

        # Input Capacitance
        try:
            z = para['Elec']['tab']['Cp'].to_numpy('float').flatten(order='F')
            zi = griddata((y, x), z, (xg, yg), method='linear')
        except:
            zi = np.zeros((len(xi), len(yi)))
        para['Elec']['tab']['Cp_2d'] = RegularGridInterpolator((xi, yi), zi, bounds_error=False, fill_value=None)

        # Output Capacitance
        try:
            z = para['Elec']['tab']['Cs'].to_numpy('float').flatten(order='F')
            zi = griddata((y, x), z, (xg, yg), method='linear')
        except:
            zi = np.zeros((len(xi), len(yi)))
        para['Elec']['tab']['Cs_2d'] = RegularGridInterpolator((xi, yi), zi, bounds_error=False, fill_value=None)

        # Coupling Capacitance
        try:
            z = para['Elec']['tab']['Cps'].to_numpy('float').flatten(order='F')
            zi = griddata((y, x), z, (xg, yg), method='linear')
        except:
            zi = np.zeros((len(xi), len(yi)))
        para['Elec']['tab']['Cps_2d'] = RegularGridInterpolator((xi, yi), zi, bounds_error=False, fill_value=None)

        # ------------------------------------------
        # Electrical (Magnetic Flux Dependency)
        # ------------------------------------------
        # Init
        x = (para['Elec']['vec']['Tj'].to_numpy() * np.ones(
            (len(para['Elec']['vec']['B'].to_numpy()), len(para['Elec']['vec']['Tj'].to_numpy())))).flatten(order='F')
        y = np.tile(para['Elec']['vec']['B'].to_numpy(), len(para['Elec']['vec']['Tj'].to_numpy()))
        xi = np.linspace(np.min(x), np.max(x), nInt)
        yi = np.linspace(np.min(y), np.max(y), nInt)
        xg, yg = np.meshgrid(yi, xi)

        # Mutual Inductance
        try:
            z = para['Elec']['tab']['Lm'].to_numpy('float').flatten(order='F')
            zi = griddata((y, x), z, (xg, yg), method='linear')
        except:
            zi = np.zeros((len(xi), len(yi)))
        para['Elec']['tab']['Lm_2d'] = RegularGridInterpolator((xi, yi), zi, bounds_error=False, fill_value=None)

        # Core Resistance
        try:
            z = para['Elec']['tab']['Rc1'].to_numpy('float').flatten(order='F')
            zi = griddata((y, x), z, (xg, yg), method='linear')
        except:
            zi = np.zeros((len(xi), len(yi)))
        para['Elec']['tab']['Rc1_2d'] = RegularGridInterpolator((xi, yi), zi, bounds_error=False, fill_value=None)

        # ------------------------------------------
        # Electrical (Voltage Dependency)
        # ------------------------------------------
        # Init
        x = (para['Elec']['vec']['Tj'].to_numpy() * np.ones(
            (len(para['Elec']['vec']['Vp'].to_numpy()), len(para['Elec']['vec']['Tj'].to_numpy())))).flatten(order='F')
        y = np.tile(para['Elec']['vec']['Vp'].to_numpy(), len(para['Elec']['vec']['Tj'].to_numpy()))
        xi = np.linspace(np.min(x), np.max(x), nInt)
        yi = np.linspace(np.min(y), np.max(y), nInt)
        xg, yg = np.meshgrid(yi, xi)

        # Mutual Inductance
        try:
            z = para['Elec']['tab']['Rc2'].to_numpy('float').flatten(order='F')
            zi = griddata((y, x), z, (xg, yg), method='linear')
        except:
            zi = np.zeros((len(xi), len(yi)))
        para['Elec']['tab']['Rc2_2d'] = RegularGridInterpolator((xi, yi), zi, bounds_error=False, fill_value=None)

    except:
        print("WARN: Two dimensional loss data could not be extracted")

    # ==============================================================================
    # Add Derived Parameters
    # ==============================================================================
    para['Elec']['con']['N'] = para['Elec']['con']['Np'] / para['Elec']['con']['Ns']
    para['Elec']['con']['Lk'] = para['Elec']['con']['Lp'] + para['Elec']['con']['Ls'] * para['Elec']['con']['N'] ** 2

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return para
