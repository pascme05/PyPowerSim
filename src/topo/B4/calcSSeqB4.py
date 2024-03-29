#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         calcSSeqB4
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
from src.general.helpFnc import deadTime, cbInter, con2dis
from src.pwm.oppPWM import oppPWM

# ==============================================================================
# External
# ==============================================================================
import numpy as np
from scipy import signal


#######################################################################################################################
# Function
#######################################################################################################################
def calcSSeqB4_CB(ref, t, Mi, setupPara, setupTopo):
    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Init
    # ==============================================================================
    x = {}
    c = {}
    s = {}
    xs = {}
    xsh = {}
    
    # ==============================================================================
    # Parameters
    # ==============================================================================
    fs = setupPara['PWM']['fs']
    Ts = 1/fs
    fel = setupTopo['fel']
    id = ['A', 'B']
    
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
    x['A'] = Mi * ref['A'] / np.max(np.abs(ref['A']))
    x['B'] = Mi * ref['B'] / np.max(np.abs(ref['B']))
    
    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    # ==============================================================================
    # Carrier
    # ==============================================================================
    # ------------------------------------------
    # Interleaved
    # ------------------------------------------
    if setupPara['PWM']['tri'] == "RE":
        c['A'] = signal.sawtooth(2*np.pi*fs*t, 1) * (-1)
        c['B'] = signal.sawtooth(2*np.pi*fs*(t - 0.5/fs), 1)
    elif setupPara['PWM']['tri'] == "FE":
        c['A'] = signal.sawtooth(2*np.pi*fs*(t - 0.5/fs), 0) * (-1)
        c['B'] = signal.sawtooth(2*np.pi*fs*t, 0)
    elif setupPara['PWM']['tri'] == "AM":
        c['A'] = signal.sawtooth(2*np.pi*fs*t, 1/3) * (-1)
        c['B'] = signal.sawtooth(2*np.pi*fs*(t - 0.5/fs), 1/3)
    else:
        c['A'] = signal.sawtooth(2*np.pi*fs*t, 0.5) * (-1)
        c['B'] = signal.sawtooth(2*np.pi*fs*(t - 0.5/fs), 0.5) * (-1)
    c['A'] = (2 * (c['A'] - min(c['A']))/(max(c['A'])-min(c['A']))) - 1
    c['B'] = (2 * (c['B'] - min(c['B']))/(max(c['B'])-min(c['B']))) - 1
    
    # ------------------------------------------
    # Non-Interleaved
    # ------------------------------------------
    if setupPara['PWM']['int'] == 0:
        c['B'] = c['A']

    # ==============================================================================
    # Sampling
    # ==============================================================================
    for i in range(0, len(id)):
        if setupPara['PWM']['samp'] == "RS":
            if setupPara['PWM']['upd'] == "SE":
                xs[id[i]] = con2dis(x[id[i]], t, Ts)
                xsh[id[i]] = np.roll(x[id[i]], int(len(xs[id[i]])*fel/fs))
            else:
                xs[id[i]] = con2dis(x[id[i]], t, Ts/2)
                xsh[id[i]] = np.roll(x[id[i]], int(len(xs[id[i]])*fel/fs/2))
        else:
            xs[id[i]] = x[id[i]]
            xsh[id[i]] = x[id[i]]

    # ==============================================================================
    # Intersections
    # ==============================================================================
    s['A'] = cbInter(xs['A'], c['A'], Mi, tmin)
    s['B'] = cbInter(xs['B'], c['B'], Mi, tmin)  
            
    ###################################################################################################################
    # Post-Processing
    ###################################################################################################################
    # ==============================================================================
    # Minimum Pulse Width
    # ==============================================================================
    if setupPara['PWM']['tmin'] > 0:
        for i in range(0, len(id)):
            hold = tmin
            for ii in range(0, len(s)):
                if hold >= tmin:
                    if Mi != 0:
                        s[id[i]][ii] = s[id[i]][ii]
                    hold = 0
                else:
                    s[id[i]][i] = s[id[i]][ii - 1]
                    hold = hold + 1

    # ==============================================================================
    # Dead-time
    # ==============================================================================
    if setupPara['PWM']['td'] > 0:
        s['A'] = deadTime(s['A'], td)
        s['B'] = deadTime(s['B'], td)

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return [xs, xsh, s, c]


#######################################################################################################################
# Function
#######################################################################################################################
def calcSSeqB4_FF(ref, t, Mi, _, setupTopo):
    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Parameters
    # ==============================================================================
    fel = setupTopo['fel']
    alpha = np.arccos(Mi*np.pi/4)
    
    # ==============================================================================
    # Variables
    # ==============================================================================
    c = {}
    s = {}
    x = {}
    
    ###################################################################################################################
    # Pre-Processing
    ###################################################################################################################
    # ==============================================================================
    # Carrier
    # ==============================================================================
    c['A'] = np.zeros(np.size(t))
    c['B'] = np.zeros(np.size(t))
    
    # ==============================================================================
    # Scale reference
    # ==============================================================================
    x['A'] = Mi * ref['A'] / np.max(np.abs(ref['A']))
    x['B'] = Mi * ref['B'] / np.max(np.abs(ref['B']))
    
    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    s['A'] = signal.square(2 * np.pi * fel * t + np.pi/2 - alpha, duty=0.5)
    s['B'] = signal.square(2 * np.pi * fel * t + np.pi/2 + alpha, duty=0.5)
            
    ###################################################################################################################
    # Post-Processing
    ###################################################################################################################
    # ==============================================================================
    # Output
    # ==============================================================================
    xs = x
    xsh = x
    
    ###################################################################################################################
    # Return
    ###################################################################################################################
    return [xs, xsh, s, c]


#######################################################################################################################
# Function
#######################################################################################################################
def calcSSeqB4_OPP(ref, t, Mi, setupPara, setupTopo):
    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Init
    # ==============================================================================
    x = {}
    c = {}
    s = {}
    xs = {}
    xsh = {}

    # ==============================================================================
    # Parameters
    # ==============================================================================
    fs = setupPara['PWM']['fs']
    Ts = 1 / fs
    fel = setupTopo['fel']
    Tel = 1 / fel
    q = int(fs / fel)
    kmax = 10 * q
    N = int((t[-1] - t[0]) / Tel)
    tmin = int(setupPara['PWM']['tmin'] / (t[1] - t[0]))
    td = int(setupPara['PWM']['td'] / (t[1] - t[0]))
    id = ['A', 'B']

    # ==============================================================================
    # Variables
    # ==============================================================================
    c['A'] = np.zeros(np.size(t))
    s['A'] = np.ones(np.size(t))
    ss = np.zeros(np.size(t))
    ang = np.linspace(0, np.pi * 2 * N, np.size(t))
    ang_total = []
    val_total = []

    ###################################################################################################################
    # Pre-Processing
    ###################################################################################################################
    # ==============================================================================
    # Scale reference
    # ==============================================================================
    x['A'] = Mi * ref['A'] / np.max(np.abs(ref['A']))
    x['B'] = Mi * ref['B'] / np.max(np.abs(ref['B']))

    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    # ==============================================================================
    # Fundamental Angles (0, 2pi)
    # ==============================================================================
    [ang_fun, val_fun, _] = oppPWM(kmax, q*2, Mi/4*np.pi, 4, setupTopo)

    # ==============================================================================
    # Complete Angles
    # ==============================================================================
    for i in range(0, N):
        if i == 0:
            ang_total = ang_fun
            val_total = val_fun
        else:
            ang_total = np.concatenate((ang_total, ang_fun + 2 * np.pi * i), axis=0)
            val_total = np.concatenate((val_total, val_fun), axis=0)

    # ==============================================================================
    # Switching times
    # ==============================================================================
    # ------------------------------------------
    # Switching Edges
    # ------------------------------------------
    for i in range(0, len(ang_total)):
        idx = np.argmin(abs(ang - ang_total[i]))
        c['A'][idx] = (-1) ** i

    # ------------------------------------------
    # Switching Function
    # ------------------------------------------
    ang2 = ang
    for i in range(1, len(ss)):
        if ang2[i] < np.pi:
            if c['A'][i] == 0:
                ss[i] = ss[i - 1]
            else:
                if c['A'][i] > 0:
                    ss[i] = 1
                else:
                    ss[i] = 0
        else:
            if c['A'][i] == 0:
                ss[i] = ss[i - 1]
            else:
                if c['A'][i] > 0:
                    ss[i] = -1
                else:
                    ss[i] = 0

        if ang2[i] > 2*np.pi:
            ang2 = ang2 - 2*np.pi

    # ==============================================================================
    # Two Phases
    # ==============================================================================
    # ------------------------------------------
    # Direction
    # ------------------------------------------
    for i in range(1, len(ss)):
        if ss[i] >= 0:
            s['A'][i] = 1
        elif ss[i] < 0:
            s['A'][i] = -1

    # ------------------------------------------
    # Interleaved
    # ------------------------------------------
    s['B'] = np.roll(s['A'], -int(np.floor(180 / 360 / N * len(ss))))
    c['B'] = np.roll(c['A'], -int(np.floor(180 / 360 / N * len(ss))))

    # ------------------------------------------
    # Non-Interleaved
    # ------------------------------------------
    if setupPara['PWM']['int'] == 0:
        c['B'] = c['A']

    # ==============================================================================
    # Sampling
    # ==============================================================================
    for i in range(0, len(id)):
        if setupPara['PWM']['samp'] == "RS":
            if setupPara['PWM']['upd'] == "SE":
                xs[id[i]] = con2dis(x[id[i]], t, Ts)
                xsh[id[i]] = np.roll(x[id[i]], int(len(xs[id[i]]) * fel / fs))
            else:
                xs[id[i]] = con2dis(x[id[i]], t, Ts / 2)
                xsh[id[i]] = np.roll(x[id[i]], int(len(xs[id[i]]) * fel / fs / 2))
        else:
            xs[id[i]] = x[id[i]]
            xsh[id[i]] = x[id[i]]

    ###################################################################################################################
    # Post-Processing
    ###################################################################################################################
    # ==============================================================================
    # Minimum Pulse Width
    # ==============================================================================
    if setupPara['PWM']['tmin'] > 0:
        for i in range(0, len(id)):
            hold = tmin
            for ii in range(0, len(s)):
                if hold >= tmin:
                    if Mi != 0:
                        s[id[i]][ii] = s[id[i]][ii]
                    hold = 0
                else:
                    s[id[i]][i] = s[id[i]][ii - 1]
                    hold = hold + 1

    # ==============================================================================
    # Dead-time
    # ==============================================================================
    if setupPara['PWM']['td'] > 0:
        s['A'] = deadTime(s['A'], td)
        s['B'] = deadTime(s['B'], td)

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return [xs, xsh, s, c]
