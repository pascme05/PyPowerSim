#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         calcLossCap
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
from src.general.helpFnc import zoh_easy

# ==============================================================================
# External
# ==============================================================================
import numpy as np
import pandas as pd
from scipy import interpolate, integrate

#######################################################################################################################
# Function
#######################################################################################################################
def calcLossSwi(i_G, i_T, i_D, v_T, v_D, t_Tj, para, setupPara, setupTopo):
    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Parameters
    # ==============================================================================
    fs = setupPara['PWM']['fs']

    # ==============================================================================
    # Variables
    # ==============================================================================
    Eon = np.zeros(np.size(i_T))
    Eoff = np.zeros(np.size(i_T))
    Erec = np.zeros(np.size(i_T))
    V_int = np.linspace(0,int(np.max(para['Swi']['Elec']['tab']['Vf'].values)), 100)
    Ciss_int, Eiss_2d = np.zeros((len(para['Swi']['Elec']['tab']['Tj']),len(V_int)))
    Coss_int, Eoss_2d = np.zeros((len(para['Swi']['Elec']['tab']['Tj']),len(V_int)))
    Crss_int, Erss_2d = np.zeros((len(para['Swi']['Elec']['tab']['Tj']),len(V_int)))
    
    # ==============================================================================
    # Output
    # ==============================================================================
    out = pd.DataFrame(columns=['p_T_c','p_T_s','p_T','p_D_c','p_D_s','p_D','p_L'])

    ###################################################################################################################
    # Pre-Processing
    ###################################################################################################################
    # ==============================================================================
    # Extract Parameters
    # ==============================================================================
    # ------------------------------------------
    # Constant
    # ------------------------------------------
    if setupPara['Elec']['SwiMdl'] == "con":
        # IGBT
        if setupPara['Elec']['SwiType'] == "IGBT":  
            Eon  = para['Swi']['Elec']['con']['Eon'] * np.ones(np.size(i_T))
            Eoff = para['Swi']['Elec']['con']['Eoff'] * np.ones(np.size(i_T))
            Erec = para['Swi']['Elec']['con']['Erec'] * np.ones(np.size(i_T))

        # MOSFET
        if setupPara['Elec']['SwiType'] == "MOSFET": 
            Eon  = 0.5 * para['Swi']['Elec']['con']['Ciss'] * np.max(np.abs(v_T))**2 * np.ones(np.size(i_T))
            Eoff = 0.5 * para['Swi']['Elec']['con']['Coss'] * np.max(np.abs(v_T))**2 * np.ones(np.size(i_T))
            Erec = 0.5 * para['Swi']['Elec']['con']['Crss'] * np.max(np.abs(v_T))**2 * np.ones(np.size(i_T))

    # ------------------------------------------
    # Piece-wise-linear
    # ------------------------------------------
    if setupPara['Elec']['SwiMdl'] == "pwl":
        # IGBT
        if setupPara['Elec']['SwiType'] == "IGBT":  
            Eon  = para['Swi']['Elec']['con']['Eon'] * i_T / para['Swi']['Elec']['con']['Imax'] * np.max(v_T) / para['Swi']['Elec']['con']['Vmax']
            Eoff = para['Swi']['Elec']['con']['Eoff'] * i_T / para['Swi']['Elec']['con']['Imax'] * np.max(v_T) / para['Swi']['Elec']['con']['Vmax']
            Erec = para['Swi']['Elec']['con']['Erec'] * i_T / para['Swi']['Elec']['con']['Imax'] * np.max(v_T) / para['Swi']['Elec']['con']['Vmax']

        # MOSFET
        if setupPara['Elec']['SwiType'] == "MOSFET": 
            Eon  = 0.5 * np.max(np.abs(v_T)) * i_T * (para['Swi']['Elec']['con']['tr'] + para['Swi']['Elec']['con']['tf']) + para['Swi']['Elec']['con']['Qrr'] * np.max(np.abs(v_T)) 
            Eoff = 0.5 * np.max(np.abs(v_T)) * i_T * (para['Swi']['Elec']['con']['tr'] + para['Swi']['Elec']['con']['tf'])
            Erec = 0.25 * para['Swi']['Elec']['con']['Qrr'] * np.max(np.abs(v_D)) * np.ones(np.size(i_T)) 

    # ------------------------------------------
    # Tabular
    # ------------------------------------------
    if setupPara['Elec']['SwiMdl'] == "tab":
        # IGBT
        if setupPara['Elec']['SwiType'] == "IGBT":   
            Eon_2d  = interpolate.interp2d(para['Swi']['Elec']['tab']['Tj'].to_numpy(), para['Swi']['Elec']['tab']['If'].to_numpy(), para['Swi']['Elec']['tab']['Eon'].to_numpy(), kind='linear')
            Eoff_2d = interpolate.interp2d(para['Swi']['Elec']['tab']['Tj'].to_numpy(), para['Swi']['Elec']['tab']['If'].to_numpy(), para['Swi']['Elec']['tab']['Eoff'].to_numpy(), kind='linear')
            Erec_2d = interpolate.interp2d(para['Swi']['Elec']['tab']['Tj'].to_numpy(), para['Swi']['Elec']['tab']['If'].to_numpy(), para['Swi']['Elec']['tab']['Erec'].to_numpy(), kind='linear')

            for i in range(0, len(i_T)):
                Eon[i]  = Eon_2d(t_Tj, i_T[i])
                Eoff[i] = Eoff_2d(t_Tj, i_T[i])
                Erec[i] = Erec_2d(t_Tj, i_D[i])

        # MOSFET
        if setupPara['Elec']['SwiType'] == "MOSFET":
            Ciss_2d = interpolate.interp1d(para['Swi']['Elec']['tab']['Tj'].to_numpy(),para['Swi']['Elec']['tab']['Vf'].to_numpy(),para['Swi']['Elec']['tab']['Ciss'].to_numpy()) 
            Coss_2d = interpolate.interp1d(para['Swi']['Elec']['tab']['Tj'].to_numpy(),para['Swi']['Elec']['tab']['Vf'].to_numpy(),para['Swi']['Elec']['tab']['Coss'].to_numpy()) 
            Crss_2d = interpolate.interp1d(para['Swi']['Elec']['tab']['Tj'].to_numpy(),para['Swi']['Elec']['tab']['Vf'].to_numpy(),para['Swi']['Elec']['tab']['Crss'].to_numpy()) 

            for i in range(0, len(para['Swi']['Elec']['tab']['Tj'])):
                for ii in range(0, len(V_int)):
                    Ciss_int[i,ii] = Ciss_2d(para['Swi']['Elec']['tab']['Tj'][i],V_int[ii])
                    Coss_int[i,ii] = Coss_2d(para['Swi']['Elec']['tab']['Tj'][i],V_int[ii])
                    Crss_int[i,ii] = Crss_2d(para['Swi']['Elec']['tab']['Tj'][i],V_int[ii])
                Qiss = integrate.cumulative_trapezoid(Ciss_int[i,:], x=V_int)
                Qoss = integrate.cumulative_trapezoid(Coss_int[i,:], x=V_int)
                Qrss = integrate.cumulative_trapezoid(Crss_int[i,:], x=V_int)
                Eiss_2d[i,:] = integrate.cumulative_trapezoid(Qiss[i,:], x=V_int)
                Eoss_2d[i,:] = integrate.cumulative_trapezoid(Qoss[i,:], x=V_int)
                Erss_2d[i,:] = integrate.cumulative_trapezoid(Qrss[i,:], x=V_int)

            for i in range(0, len(i_T)):
                Eon[i]  = Eiss_2d(t_Tj, i_T[i])
                Eoff[i] = Eiss_2d(t_Tj, i_T[i])
                Erec[i] = Eiss_2d(t_Tj, i_D[i])

    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    # ==============================================================================
    # Transistor
    # ==============================================================================
    # ------------------------------------------
    # Conduction
    # ------------------------------------------
    out['p_T_c'] = i_T * v_T

    # ------------------------------------------
    # Switching
    # ------------------------------------------
    if setupPara['Elec']['SwiMdl'] == "tab" or setupPara['Elec']['SwiMdl'] == "pwl":
        out['p_T_s'] = (zoh_easy(Eon, np.diff(i_G,prepend=0)) + np.roll(zoh_easy(Eoff, np.diff(i_G,append=0)*(-1)),1))*fs
    else:
        out['p_T_s'] = (Eon + Eoff)*fs

    # ==============================================================================
    # Diode
    # ==============================================================================
    # ------------------------------------------
    # Conduction
    # ------------------------------------------
    out['p_D_c'] = i_D * v_D
    
    # ------------------------------------------
    # Switching
    # ------------------------------------------
    if setupPara['Elec']['SwiMdl'] == "tab" or setupPara['Elec']['SwiMdl'] == "pwl":
        out['p_D_s'] = zoh_easy(Erec, np.diff(i_G,append=0)*(-1)) * fs
    else:
        out['p_D_s'] = Erec*fs

    # ==============================================================================
    # Total
    # ==============================================================================
    out['p_T'] = out['p_T_c'] + out['p_T_s']
    out['p_D'] = out['p_D_c'] + out['p_D_s']
    out['p_L'] = out['p_T'] + out['p_D']

    ###################################################################################################################
    # Post-Processing
    ###################################################################################################################
    # ==============================================================================
    # Scale Number Switches
    # ==============================================================================
    out['p_L'] = out['p_L'] * setupPara['Elec']['SwiSeries'] * setupPara['Elec']['SwiPara']

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return out