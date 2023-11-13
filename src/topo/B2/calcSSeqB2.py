#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         calcSSeqB2
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
from src.general.helpFnc import cbInter, con2dis, deadTime

# ==============================================================================
# External
# ==============================================================================
import numpy as np
from scipy import signal


#######################################################################################################################
# Function
#######################################################################################################################
def calcSSeqB2_CB(ref, t, Mi, setupPara, setupTopo):
    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Parameters
    # ==============================================================================
    fs = setupPara['PWM']['fs']
    Ts = 1/fs
    fel = setupTopo['fel']
    
    # ==============================================================================
    # Variables
    # ==============================================================================
    tmin = int(setupPara['PWM']['tmin']/(t[1]-t[0]))
    td = int(setupPara['PWM']['td']/(t[1]-t[0]))

    ###################################################################################################################
    # Pre-Processing
    ###################################################################################################################
    # ==============================================================================
    # Scale reference
    # ==============================================================================
    x = Mi * ref / np.max(ref)
    
    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    # ==============================================================================
    # Carrier
    # ==============================================================================
    if setupPara['PWM']['tri'] == "RE":
        c = signal.sawtooth(2*np.pi*fs*t, 1) * (-1)
    elif setupPara['PWM']['tri'] == "FE":
        c = signal.sawtooth(2*np.pi*fs*(t - 0.5/fs), 0) * (-1)
    elif setupPara['PWM']['tri'] == "AM":
        c = signal.sawtooth(2*np.pi*fs*t, 1/3) * (-1)
    else:
        c = signal.sawtooth(2*np.pi*fs*t, 0.5) * (-1)
    c = (2 * (c - min(c))/(max(c)-min(c))) - 1

    # ==============================================================================
    # Sampling
    # ==============================================================================
    if setupPara['PWM']['samp'] == "RS":
        if setupPara['PWM']['upd'] == "SE":
            xs = con2dis(x, t, Ts)
            xsh = np.roll(x, int(len(xs)*fel/fs))
        else:
            xs = con2dis(x, t, Ts/2)
            xsh = np.roll(x, int(len(xs)*fel/fs/2))
    else:
        xs = x
        xsh = x

    # ==============================================================================
    # Intersections
    # ==============================================================================
    s = cbInter(xs, c, Mi, tmin)
            
    ###################################################################################################################
    # Post-Processing
    ###################################################################################################################
    # ==============================================================================
    # Dead-time
    # ==============================================================================
    if setupPara['PWM']['td'] > 0:
        s = deadTime(s, td)

    # ==============================================================================
    # Minimum Pulse Width
    # ==============================================================================
    if setupPara['PWM']['tmin'] > 0:
        hold = tmin
        for i in range(0, len(s)):
            if hold >= tmin:
                if Mi != 0:
                    s[i] = s[i]
                hold = 0
            else:
                s[i] = s[i - 1]
                hold = hold + 1
    
    ###################################################################################################################
    # Return
    ###################################################################################################################
    return [xs, xsh, s, c]


#######################################################################################################################
# Function
#######################################################################################################################
def calcSSeqB2_FF(ref, t, Mi, _, setupTopo):
    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Init
    # ==============================================================================
    
    # ==============================================================================
    # Parameters
    # ==============================================================================
    fel = setupTopo['fel']
    
    ###################################################################################################################
    # Pre-Processing
    ###################################################################################################################
    # ==============================================================================
    # Carrier
    # ==============================================================================
    c = np.zeros(np.size(t))
    
    # ==============================================================================
    # Scale reference
    # ==============================================================================
    x = Mi * ref / np.max(ref)
    
    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    s = signal.square(2 * np.pi * fel * t, duty=0.5)
            
    ###################################################################################################################
    # Post-Processing
    ###################################################################################################################
    xs = x
    xsh = x
    
    ###################################################################################################################
    # Return
    ###################################################################################################################
    return [xs, xsh, s, c]
