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
def calcLossSwi(i_G, i_T, i_D, v_T, v_D, t_Tj, para, setupPara, setupExp):
    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Parameters
    # ==============================================================================
    fs = setupPara['PWM']['fs']
    nInt = setupExp['int']

    # ==============================================================================
    # Variables
    # ==============================================================================
    Eon = np.zeros(np.size(i_T))
    Eoff = np.zeros(np.size(i_T))
    Erec = np.zeros(np.size(i_T))
    tf = np.zeros(np.size(i_T))
    tr = np.zeros(np.size(i_T))
    try:
        V_int = np.linspace(0,int(np.max(para['Swi']['Elec']['vec']['Vf'].values)), nInt)
    except:
        V_int = np.linspace(0,int(np.max(np.abs(v_D))), nInt)
    Coss_int = np.zeros((len(para['Swi']['Elec']['vec']['Tj']),len(V_int)))
    Crss_int = np.zeros((len(para['Swi']['Elec']['vec']['Tj']),len(V_int)))
    Eoss_2d = np.zeros((len(para['Swi']['Elec']['vec']['Tj']),len(V_int)))
    Erss_2d = np.zeros((len(para['Swi']['Elec']['vec']['Tj']),len(V_int)))
    
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
            Eon  = para['Swi']['Elec']['con']['Eon']  * np.ones(np.size(i_T)) * np.max(np.abs(i_T)) / para['Swi']['Elec']['con']['Imax'] * np.max(np.abs(v_T)) / para['Swi']['Elec']['con']['Vmax']
            Eoff = para['Swi']['Elec']['con']['Eoff'] * np.ones(np.size(i_T)) * np.max(np.abs(i_T)) / para['Swi']['Elec']['con']['Imax'] * np.max(np.abs(v_T)) / para['Swi']['Elec']['con']['Vmax']
            Erec = para['Swi']['Elec']['con']['Erec'] * np.ones(np.size(i_T)) * np.max(np.abs(i_T)) / para['Swi']['Elec']['con']['Imax'] * np.max(np.abs(v_T)) / para['Swi']['Elec']['con']['Vmax']

        # MOSFET
        if setupPara['Elec']['SwiType'] == "MOSFET": 
            Eon  = 0.50 * np.max(np.abs(v_T)) * np.max(np.abs(i_T)) * (para['Swi']['Elec']['con']['tr'] + para['Swi']['Elec']['con']['tf']) + 0.5 * para['Swi']['Elec']['con']['Coss'] * np.max(np.abs(v_T))**2 * np.ones(np.size(i_T))
            Eoff = 0.50 * np.max(np.abs(v_T)) * np.max(np.abs(i_T)) * (para['Swi']['Elec']['con']['tr'] + para['Swi']['Elec']['con']['tf']) + 0.5 * para['Swi']['Elec']['con']['Coss'] * np.max(np.abs(v_T))**2 * np.ones(np.size(i_T))
            Erec = 0.25 * para['Swi']['Elec']['con']['Qrr'] * np.max(np.abs(v_D)) * np.ones(np.size(i_T))  + 0.5 * para['Swi']['Elec']['con']['Crss'] * np.max(np.abs(v_T))**2 * np.ones(np.size(i_T))

    # ------------------------------------------
    # Piece-wise-linear
    # ------------------------------------------
    if setupPara['Elec']['SwiMdl'] == "pwl":
        # IGBT
        if setupPara['Elec']['SwiType'] == "IGBT":  
            Eon  = para['Swi']['Elec']['con']['Eon']  * np.abs(i_T) / para['Swi']['Elec']['con']['Imax'] * np.max(np.abs(v_T)) / para['Swi']['Elec']['con']['Vmax']
            Eoff = para['Swi']['Elec']['con']['Eoff'] * np.abs(i_T) / para['Swi']['Elec']['con']['Imax'] * np.max(np.abs(v_T)) / para['Swi']['Elec']['con']['Vmax']
            Erec = para['Swi']['Elec']['con']['Erec'] * np.abs(i_T) / para['Swi']['Elec']['con']['Imax'] * np.max(np.abs(v_T)) / para['Swi']['Elec']['con']['Vmax']

        # MOSFET
        if setupPara['Elec']['SwiType'] == "MOSFET":
            Eon  = 0.50 * np.max(np.abs(v_T)) * np.abs(i_T) * (para['Swi']['Elec']['con']['tr'] + para['Swi']['Elec']['con']['tf']) + para['Swi']['Elec']['con']['Qrr'] * np.max(np.abs(v_T)) 
            Eoff = 0.50 * np.max(np.abs(v_T)) * np.abs(i_T) * (para['Swi']['Elec']['con']['tr'] + para['Swi']['Elec']['con']['tf'])
            Erec = 0.25 * para['Swi']['Elec']['con']['Qrr'] * np.max(np.abs(v_D)) * np.ones(np.size(i_T)) 

    # ------------------------------------------
    # Tabular
    # ------------------------------------------
    if setupPara['Elec']['SwiMdl'] == "tab":
        # IGBT
        if setupPara['Elec']['SwiType'] == "IGBT" or setupPara['PWM']['swloss'] == 0:   
            Eon_2d  = interpolate.interp2d(para['Swi']['Elec']['vec']['Tj'].to_numpy(), para['Swi']['Elec']['vec']['If'].to_numpy(), para['Swi']['Elec']['tab']['Eon'].to_numpy(), kind='linear')
            Eoff_2d = interpolate.interp2d(para['Swi']['Elec']['vec']['Tj'].to_numpy(), para['Swi']['Elec']['vec']['If'].to_numpy(), para['Swi']['Elec']['tab']['Eoff'].to_numpy(), kind='linear')
            Erec_2d = interpolate.interp2d(para['Swi']['Elec']['vec']['Tj'].to_numpy(), para['Swi']['Elec']['vec']['If'].to_numpy(), para['Swi']['Elec']['tab']['Erec'].to_numpy(), kind='linear')

            for i in range(0, len(i_T)):
                Eon[i]  = Eon_2d(t_Tj, i_T[i])
                Eoff[i] = Eoff_2d(t_Tj, i_T[i])
                Erec[i] = Erec_2d(t_Tj, i_D[i])

        # MOSFET
        if setupPara['Elec']['SwiType'] == "MOSFET" and setupPara['PWM']['swloss'] == 1:
            Coss_2d = interpolate.interp2d(para['Swi']['Elec']['vec']['Tj'].to_numpy(),para['Swi']['Elec']['vec']['Vf'].to_numpy(),para['Swi']['Elec']['tab']['Coss'].to_numpy(), kind='linear') 
            Crss_2d = interpolate.interp2d(para['Swi']['Elec']['vec']['Tj'].to_numpy(),para['Swi']['Elec']['vec']['Vf'].to_numpy(),para['Swi']['Elec']['tab']['Crss'].to_numpy(), kind='linear') 

            IGon  = (para['Swi']['Elec']['con']['Vg'] - para['Swi']['Elec']['con']['Vpl']) / para['Swi']['Elec']['con']['Rg']
            IGoff = para['Swi']['Elec']['con']['Vpl'] / para['Swi']['Elec']['con']['Rg']

            for i in range(0, len(para['Swi']['Elec']['vec']['Tj'])):
                for ii in range(0, len(V_int)):
                    Coss_int[i,ii] = Coss_2d(para['Swi']['Elec']['vec']['Tj'][i],V_int[ii])
                    Crss_int[i,ii] = Crss_2d(para['Swi']['Elec']['vec']['Tj'][i],V_int[ii])
                Qoss = integrate.cumulative_trapezoid(Coss_int[i,:], x=V_int, initial=V_int[0])
                Qrss = integrate.cumulative_trapezoid(Crss_int[i,:], x=V_int, initial=V_int[0])
                Eoss_2d[i,:] = integrate.cumulative_trapezoid(Qoss, x=V_int, initial=V_int[0])
                Erss_2d[i,:] = integrate.cumulative_trapezoid(Qrss, x=V_int, initial=V_int[0])

            Eon_2d  = interpolate.interp2d(para['Swi']['Elec']['vec']['Tj'].to_numpy(), V_int, np.transpose(Eoss_2d), kind='linear')
            Eoff_2d = interpolate.interp2d(para['Swi']['Elec']['vec']['Tj'].to_numpy(), V_int, np.transpose(Eoss_2d), kind='linear')
            Erec_2d = interpolate.interp2d(para['Swi']['Elec']['vec']['Tj'].to_numpy(), V_int, np.transpose(Erss_2d), kind='linear')

            for i in range(0, len(i_T)):
                Eon[i]  = Eon_2d(t_Tj, abs(v_T[i]))
                Eoff[i] = Eoff_2d(t_Tj, abs(v_T[i]))
                Erec[i] = Erec_2d(t_Tj, abs(v_D[i]))
                tf1 = (np.max(np.abs(v_T)) - np.abs(v_T[i]))*(Crss_2d(t_Tj,np.abs(v_T[i]))/IGon)
                tf2 = (np.max(np.abs(v_T)) - np.abs(v_T[i]))*(Crss_2d(t_Tj,np.max(np.abs(v_T)))/IGon)
                tr1 = (np.max(np.abs(v_T)) - np.abs(v_T[i]))*(Crss_2d(t_Tj,np.abs(v_T[i]))/IGoff)
                tr2 = (np.max(np.abs(v_T)) - np.abs(v_T[i]))*(Crss_2d(t_Tj,np.max(np.abs(v_T)))/IGoff)
                tf[i] = (tf1 + tf2)/2
                tr[i] = (tr1 + tr2)/2

            Eon = Eon + 0.5*np.abs(i_T)*np.max(abs(v_T))*(tr + tf)
            Eoff = Eoff + 0.5*np.abs(i_T)*np.max(abs(v_T))*(tr + tf)
            Erec = Erec + 0.25 * para['Swi']['Elec']['con']['Qrr'] * np.max(np.abs(v_D))

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
    out['p_T_s'] = (zoh_easy(Eon, np.diff(i_G,prepend=0)) + np.roll(zoh_easy(Eoff, np.diff(i_G,append=0)*(-1)),1))*fs
    ## IGBT
    #if setupPara['Elec']['SwiType'] == "IGBT":
    #    out['p_T_s'] = (zoh_easy(Eon, np.diff(i_G,prepend=0)) + np.roll(zoh_easy(Eoff, np.diff(i_G,append=0)*(-1)),1))*fs

    ## MOSFET
    #if setupPara['Elec']['SwiType'] == "MOSFET":
     #   if setupPara['Elec']['SwiMdl'] == "tab":
      #      if setupPara['PWM']['swloss'] == 0:
      #          out['p_T_s'] = (zoh_easy(Eon, np.diff(i_G,prepend=0)) + np.roll(zoh_easy(Eoff, np.diff(i_G,append=0)*(-1)),1))*fs
      #      else:
      #          out['p_T_s'] = 0.5*np.abs(i_T)*np.max(abs(v_T))*(tr + tf)*fs + (zoh_easy(Eon, np.diff(i_G,prepend=0)) + np.roll(zoh_easy(Eoff, np.diff(i_G,append=0)*(-1)),1))*fs
      #  elif setupPara['Elec']['SwiMdl'] == "pwl":
       #     out['p_T_s'] = (zoh_easy(Eon, np.diff(i_G,prepend=0)) + np.roll(zoh_easy(Eoff, np.diff(i_G,append=0)*(-1)),1))*fs
      #  else:
            #out['p_T_s'] = (Eon + Eoff)*fs
       #     out['p_T_s'] = (zoh_easy(Eon, np.diff(i_G,prepend=0)) + np.roll(zoh_easy(Eoff, np.diff(i_G,append=0)*(-1)),1))*fs

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
    out['p_D_s'] = zoh_easy(Erec, np.diff(i_G,append=0)*(-1)) * fs

    # IGBT
    #if setupPara['Elec']['SwiType'] == "IGBT":
    #    out['p_D_s'] = zoh_easy(Erec, np.diff(i_G,append=0)*(-1)) * fs

    # MOSFET
    #if setupPara['Elec']['SwiType'] == "MOSFET":
     #   if setupPara['Elec']['SwiMdl'] == "tab":
       #     if setupPara['PWM']['swloss'] == 0:
       #         out['p_D_s'] = zoh_easy(Erec, np.diff(i_G,append=0)*(-1)) * fs
       #     else:
       #         out['p_D_s'] = zoh_easy(Erec, np.diff(i_G,append=0)*(-1)) * fs + 0.25 * para['Swi']['Elec']['con']['Qrr'] * np.max(np.abs(v_D)) * np.ones(np.size(i_T))
       # elif setupPara['Elec']['SwiMdl'] == "pwl":
      #      out['p_D_s'] = zoh_easy(Erec, np.diff(i_G,append=0)*(-1)) * fs 
       # else:
            #out['p_D_s'] = Erec*fs
        #    out['p_D_s'] = zoh_easy(Erec, np.diff(i_G,append=0)*(-1)) * fs

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