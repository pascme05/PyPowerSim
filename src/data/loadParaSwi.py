#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         loadParaSwi
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
This function loads the parameters for the switching device under \par. This includes experimental, data, topology, and
electrical as well as thermal parameter information. The parameters are summarized in one common para variable.
Inputs:     1) name:    name of the parameter file for the switches
            2) path:    includes all path variables
            2) setup:   includes all simulation variables
Outputs:    1) para:    output parameter file for the switches
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
from scipy import integrate
from scipy.interpolate import griddata, RegularGridInterpolator


#######################################################################################################################
# Function
#######################################################################################################################
def loadParaSwi(name, path, setup):
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
    lenElec = 24
    lenTher = 6

    ###################################################################################################################
    # Loading Data
    ###################################################################################################################
    # ==============================================================================
    # Path
    # ==============================================================================
    name = name + '.xlsx'
    filename = pjoin(path, 'Swi', name)

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
    para['Elec']['tab']['Vf'] = para['Elec']['tab']['Vf'].dropna(axis=0, how='all')
    para['Elec']['tab']['Vf'] = para['Elec']['tab']['Vf'].dropna(axis=1, how='all')
    para['Elec']['tab']['Eon'] = para['Elec']['tab']['Eon'].dropna(axis=0, how='all')
    para['Elec']['tab']['Eon'] = para['Elec']['tab']['Eon'].dropna(axis=1, how='all')
    para['Elec']['tab']['Eoff'] = para['Elec']['tab']['Eoff'].dropna(axis=0, how='all')
    para['Elec']['tab']['Eoff'] = para['Elec']['tab']['Eoff'].dropna(axis=1, how='all')
    para['Elec']['tab']['Erec'] = para['Elec']['tab']['Erec'].dropna(axis=0, how='all')
    para['Elec']['tab']['Erec'] = para['Elec']['tab']['Erec'].dropna(axis=1, how='all')
    para['Elec']['tab']['Vfd'] = para['Elec']['tab']['Vfd'].dropna(axis=0, how='all')
    para['Elec']['tab']['Vfd'] = para['Elec']['tab']['Vfd'].dropna(axis=1, how='all')
    para['Elec']['tab']['Ciss'] = para['Elec']['tab']['Ciss'].dropna(axis=0, how='all')
    para['Elec']['tab']['Ciss'] = para['Elec']['tab']['Ciss'].dropna(axis=1, how='all')
    para['Elec']['tab']['Coss'] = para['Elec']['tab']['Coss'].dropna(axis=0, how='all')
    para['Elec']['tab']['Coss'] = para['Elec']['tab']['Coss'].dropna(axis=1, how='all')
    para['Elec']['tab']['Crss'] = para['Elec']['tab']['Crss'].dropna(axis=0, how='all')
    para['Elec']['tab']['Crss'] = para['Elec']['tab']['Crss'].dropna(axis=1, how='all')

    # ==============================================================================
    # Losses
    # ==============================================================================
    if setup['Par']['PWM']['loss'] == 0:
        # ------------------------------------------
        # Matrix
        # ------------------------------------------
        para['Elec']['tab']['Vf'] = para['Elec']['tab']['Vf'] * 0
        para['Elec']['tab']['Eon'] = para['Elec']['tab']['Eon'] * 0
        para['Elec']['tab']['Eoff'] = para['Elec']['tab']['Eoff'] * 0
        para['Elec']['tab']['Erec'] = para['Elec']['tab']['Erec'] * 0
        para['Elec']['tab']['Vfd'] = para['Elec']['tab']['Vfd'] * 0
        para['Elec']['tab']['Ciss'] = para['Elec']['tab']['Ciss'] * 0
        para['Elec']['tab']['Coss'] = para['Elec']['tab']['Coss'] * 0
        para['Elec']['tab']['Crss'] = para['Elec']['tab']['Crss'] * 0

        # ------------------------------------------
        # Constant
        # ------------------------------------------
        para['Elec']['con']['Vf'] = para['Elec']['con']['Vf'] * 0
        para['Elec']['con']['Eon'] = para['Elec']['con']['Eon'] * 0
        para['Elec']['con']['Eoff'] = para['Elec']['con']['Eoff'] * 0
        para['Elec']['con']['Erec'] = para['Elec']['con']['Erec'] * 0
        para['Elec']['con']['Vfd'] = para['Elec']['con']['Vfd'] * 0
        para['Elec']['con']['Ron'] = para['Elec']['con']['Ron'] * 0
        para['Elec']['con']['Roff'] = para['Elec']['con']['Roff'] * 0
        para['Elec']['con']['RonD'] = para['Elec']['con']['RonD'] * 0
        para['Elec']['con']['RoffD'] = para['Elec']['con']['RoffD'] * 0
        para['Elec']['con']['Ciss'] = para['Elec']['con']['Ciss'] * 0
        para['Elec']['con']['Coss'] = para['Elec']['con']['Coss'] * 0
        para['Elec']['con']['Crss'] = para['Elec']['con']['Crss'] * 0

    # ==============================================================================
    # Interpolation Functions
    # ==============================================================================
    try:
        # ------------------------------------------
        # Init
        # ------------------------------------------
        # General
        nInt = setup['Exp']['int']

        # Voltages
        try:
            V_int = np.linspace(0, int(np.max(para['Elec']['vec']['Vf'].values)), nInt)
        except:
            V_int = np.linspace(0, int(np.max(np.abs(setup['Dat']['stat']['Vdc']))), nInt)

        # Energies
        Coss_int = np.zeros((len(para['Elec']['vec']['Tj']), len(V_int)))
        Crss_int = np.zeros((len(para['Elec']['vec']['Tj']), len(V_int)))
        Eoss_2d = np.zeros((len(para['Elec']['vec']['Tj']), len(V_int)))
        Erss_2d = np.zeros((len(para['Elec']['vec']['Tj']), len(V_int)))

        # ------------------------------------------
        # Elec
        # ------------------------------------------
        # Transistor
        x = (para['Elec']['vec']['Tj'].to_numpy() * np.ones((len(para['Elec']['vec']['If'].to_numpy()), len(para['Elec']['vec']['Tj'].to_numpy())))).flatten(order='F')
        y = np.tile(para['Elec']['vec']['If'].to_numpy(), len(para['Elec']['vec']['Tj'].to_numpy()))
        z = para['Elec']['tab']['Vf'].to_numpy('float').flatten(order='F')
        xi = np.linspace(np.min(x), np.max(x), nInt)
        yi = np.linspace(np.min(y), np.max(y), nInt)
        xg, yg = np.meshgrid(yi, xi)
        zi = griddata((y, x), z, (xg, yg), method='linear')
        para['Elec']['tab']['Vce_2d'] = RegularGridInterpolator((xi, yi), zi, bounds_error=False, fill_value=None)

        # Diode
        x = (para['Elec']['vec']['Tj'].to_numpy() * np.ones((len(para['Elec']['vec']['Ifd'].to_numpy()), len(para['Elec']['vec']['Tj'].to_numpy())))).flatten(order='F')
        y = np.tile(para['Elec']['vec']['Ifd'].to_numpy(), len(para['Elec']['vec']['Tj'].to_numpy()))
        z = para['Elec']['tab']['Vfd'].to_numpy('float').flatten(order='F')
        xi = np.linspace(np.min(x), np.max(x), nInt)
        yi = np.linspace(np.min(y), np.max(y), nInt)
        xg, yg = np.meshgrid(yi, xi)
        zi = griddata((y, x), z, (xg, yg), method='linear')
        para['Elec']['tab']['Vfd_2d'] = RegularGridInterpolator((xi, yi), zi, bounds_error=False, fill_value=None)

        # ------------------------------------------
        # Losses
        # ------------------------------------------
        # IGBT
        if setup['Par']['Elec']['SwiType'] == "IGBT" or setup['Par']['PWM']['swloss'] == 0:
            x = (para['Elec']['vec']['Tj'].to_numpy() * np.ones((len(para['Elec']['vec']['If'].to_numpy()), len(para['Elec']['vec']['Tj'].to_numpy())))).flatten(order='F')
            y = np.tile(para['Elec']['vec']['If'].to_numpy(), len(para['Elec']['vec']['Tj'].to_numpy()))
            z = para['Elec']['tab']['Eon'].to_numpy('float').flatten(order='F')
            xi = np.linspace(np.min(x), np.max(x), nInt)
            yi = np.linspace(np.min(y), np.max(y), nInt)
            xg, yg = np.meshgrid(yi, xi)
            zi = griddata((y, x), z, (xg, yg), method='linear')
            para['Elec']['tab']['Eon_2d'] = RegularGridInterpolator((xi, yi), zi, bounds_error=False, fill_value=None)

            z = para['Elec']['tab']['Eoff'].to_numpy('float').flatten(order='F')
            zi = griddata((y, x), z, (xg, yg), method='linear')
            para['Elec']['tab']['Eoff_2d'] = RegularGridInterpolator((xi, yi), zi, bounds_error=False, fill_value=None)

            z = para['Elec']['tab']['Erec'].to_numpy('float').flatten(order='F')
            zi = griddata((y, x), z, (xg, yg), method='linear')
            para['Elec']['tab']['Erec_2d'] = RegularGridInterpolator((xi, yi), zi, bounds_error=False, fill_value=None)

        # MOSFET
        if setup['Par']['Elec']['SwiType'] == "MOSFET" and setup['Par']['PWM']['swloss'] == 1:
            x = (para['Elec']['vec']['Tj'].to_numpy() * np.ones((len(para['Elec']['vec']['Vf'].to_numpy()), len(para['Elec']['vec']['Tj'].to_numpy())))).flatten(order='F')
            y = np.tile(para['Elec']['vec']['Vf'].to_numpy(), len(para['Elec']['vec']['Tj'].to_numpy()))
            z = para['Elec']['tab']['Coss'].to_numpy('float').flatten(order='F')
            xi = np.linspace(np.min(x), np.max(x), nInt)
            yi = np.linspace(np.min(y), np.max(y), nInt)
            xg, yg = np.meshgrid(yi, xi)
            zi = griddata((y, x), z, (xg, yg), method='linear')
            para['Elec']['tab']['Coss_2d'] = RegularGridInterpolator((xi, yi), zi, bounds_error=False, fill_value=None)

            z = para['Elec']['tab']['Crss'].to_numpy('float').flatten(order='F')
            zi = griddata((y, x), z, (xg, yg), method='linear')
            para['Elec']['tab']['Crss_2d'] = RegularGridInterpolator((xi, yi), zi, bounds_error=False, fill_value=None)

            for i in range(0, len(para['Elec']['vec']['Tj'])):
                for ii in range(0, len(V_int)):
                    Coss_int[i, ii] = para['Elec']['tab']['Coss_2d'](para['Elec']['vec']['Tj'][i], V_int[ii])
                    Crss_int[i, ii] = para['Elec']['tab']['Crss_2d'](para['Elec']['vec']['Tj'][i], V_int[ii])
                Qoss = integrate.cumulative_trapezoid(Coss_int[i, :], x=V_int, initial=V_int[0])
                Qrss = integrate.cumulative_trapezoid(Crss_int[i, :], x=V_int, initial=V_int[0])
                Eoss_2d[i, :] = integrate.cumulative_trapezoid(Qoss, x=V_int, initial=V_int[0])
                Erss_2d[i, :] = integrate.cumulative_trapezoid(Qrss, x=V_int, initial=V_int[0])

            x = (para['Elec']['vec']['Tj'].to_numpy() * np.ones((len(V_int), len(para['Elec']['vec']['Tj'].to_numpy())))).flatten(order='F')
            y = np.tile(V_int, len(para['Elec']['vec']['Tj'].to_numpy()))
            z = np.transpose(Eoss_2d).flatten(order='F')
            xi = np.linspace(np.min(x), np.max(x), nInt)
            yi = np.linspace(np.min(y), np.max(y), nInt)
            xg, yg = np.meshgrid(yi, xi)
            zi = griddata((y, x), z, (xg, yg), method='linear')
            para['Elec']['tab']['Eon_2d'] = RegularGridInterpolator((xi, yi), zi, bounds_error=False, fill_value=None)

            z = np.transpose(Eoss_2d).flatten(order='F')
            zi = griddata((y, x), z, (xg, yg), method='linear')
            para['Elec']['tab']['Eoff_2d'] = RegularGridInterpolator((xi, yi), zi, bounds_error=False, fill_value=None)

            z = np.transpose(Erss_2d).flatten(order='F')
            zi = griddata((y, x), z, (xg, yg), method='linear')
            para['Elec']['tab']['Erec_2d'] = RegularGridInterpolator((xi, yi), zi, bounds_error=False, fill_value=None)
    except:
        print("WARN: Two dimensional loss data could not be extracted")

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return para
