#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         calcDistB2
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
import numpy as np
from scipy.fft import fft

#######################################################################################################################
# Function
#######################################################################################################################
def calcDistB2_Num(t, i_a, v_a, i_dc, v_dc, Vdc, setupTopo):
    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Variables
    # ==============================================================================
    distAc = {}
    distDc = {}
    
    # ==============================================================================
    # Parameters
    # ==============================================================================
    fel = setupTopo['fel']
    Tel = 1/fel
    dt = t[1] - t[0]
    K = int(np.round((t[-1]-t[0])/Tel))
    N = int(len(v_a)) 
    
    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    # ==============================================================================
    # AC-Side
    # ==============================================================================
    V_a0_eff = np.sqrt(1/Tel/K * np.sum(v_a**2 * dt))
    V_a0_v1_eff = (1/np.sqrt(2))*2*np.abs(fft(v_a)/N)[K]
    V_a0_thd = np.sqrt(V_a0_eff**2 - V_a0_v1_eff**2)/V_a0_eff * Vdc/2
    I_a_eff = np.sqrt(1/Tel/K * np.sum(i_a**2 * dt))
    I_a_v1_eff = (1/np.sqrt(2))*2*np.abs(fft(i_a)/N)[K]
    I_a_thd = np.sqrt(I_a_eff**2 - I_a_v1_eff**2) 
    
    # ==============================================================================
    # DC-Side
    # ==============================================================================
    V_dc_eff = np.sqrt(1/Tel/K * np.sum(v_dc**2 * dt))
    V_dc_v1_eff = np.abs(fft(v_dc)/N)[0]
    V_dc_thd = np.sqrt(V_dc_eff**2 - V_dc_v1_eff**2)
    I_dc_eff = np.sqrt(1/Tel/K * np.sum(i_dc**2 * dt))
    I_dc_v1_eff = np.abs(fft(i_dc)/N)[0]
    I_dc_thd = np.sqrt(I_dc_eff**2 - I_dc_v1_eff**2)
        
    ###################################################################################################################
    # Post-Processing
    ###################################################################################################################
    # ==============================================================================
    # AC-Side
    # ==============================================================================
    distAc['V_a0_eff'] = V_a0_eff
    distAc['V_a0_v1_eff'] = V_a0_v1_eff
    distAc['V_a0_thd'] = V_a0_thd
    distAc['I_a_eff'] = I_a_eff
    distAc['I_a_v1_eff'] = I_a_v1_eff
    distAc['I_a_thd'] = I_a_thd
    
    # ==============================================================================
    # DC-Side
    # ==============================================================================
    distDc['V_dc_eff'] = V_dc_eff
    distDc['V_dc_v1_eff'] = V_dc_v1_eff
    distDc['V_dc_thd'] = V_dc_thd
    distDc['I_dc_eff'] = I_dc_eff
    distDc['I_dc_v1_eff'] = I_dc_v1_eff
    distDc['I_dc_thd'] = I_dc_thd
    
    ###################################################################################################################
    # Return
    ###################################################################################################################
    return [distAc, distDc]


#######################################################################################################################
# Function
#######################################################################################################################
def calcDistB2_Ana(Mi, Vdc, setupTopo, setupPara):
    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Init
    # ==============================================================================
    distAc = {}
    distDc = {}
    
    # ==============================================================================
    # Parameters
    # ==============================================================================
    fs = setupPara['PWM']['fs']
    fel = setupTopo['fel']
    Ts = 1/fs
    L = setupTopo['L']
    R = setupTopo['R']
    Z = np.sqrt(R**2 + (2*np.pi*fel*L)**2)
    
    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    # ==============================================================================
    # AC-Side
    # ==============================================================================
    V_a0_eff = Vdc/2
    V_a0_v1_eff = 1/np.sqrt(2) * Vdc/2 * Mi
    V_a0_thd = Vdc/2 * np.sqrt(1 - Mi**2/2)
    I_a_eff = V_a0_eff / Z 
    I_a_v1_eff = V_a0_v1_eff / Z 
    I_a_thd = 1/np.sqrt(48) * Vdc/2 * Ts/L * np.sqrt(3/8*Mi**4 - Mi**2 +1)
    
    # ==============================================================================
    # DC-Side
    # ==============================================================================
    V_dc_eff = Mi * 4 / np.pi * Vdc/2
    I_dc_thd = 1/np.sqrt(48) * Vdc/2 * Ts/L * np.sqrt(3/8*Mi**4 - Mi**2 +1)
    
    ###################################################################################################################
    # Post-Processing
    ###################################################################################################################
    # ==============================================================================
    # AC-Side
    # ==============================================================================
    distAc['V_a0_eff'] = V_a0_eff
    distAc['V_a0_v1_eff'] = V_a0_v1_eff
    distAc['V_a0_thd'] = V_a0_thd
    distAc['I_a_eff'] = I_a_eff
    distAc['I_a_v1_eff'] = I_a_v1_eff
    distAc['I_a_thd'] = I_a_thd
    
    # ==============================================================================
    # DC-Side
    # ==============================================================================
    distDc['V_dc_eff'] = V_dc_eff
    distDc['V_dc_v1_eff'] = 0
    distDc['V_dc_thd'] = 0
    distDc['I_dc_eff'] = 0
    distDc['I_dc_v1_eff'] = 0
    distDc['I_dc_thd'] = I_dc_thd
    
    ###################################################################################################################
    # Return
    ###################################################################################################################
    return [distAc, distDc]