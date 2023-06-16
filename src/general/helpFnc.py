#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         helpFnc
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
import math
from scipy.fft import fft
from os.path import dirname, join as pjoin
import os

#######################################################################################################################
# Init Setup files
#######################################################################################################################
def initSetup():
    setupExp = {}
    setupData = {}
    setupData['stat'] = {}
    setupData['trans'] = {}
    setupTopo = {}
    setupPara = {}
    setupPara['PWM'] = {}
    setupPara['Elec'] = {}
    setupPara['Ther'] = {}
    setupPara['Stat'] = {}
    setupPara['Life'] = {}

    return [setupExp, setupData, setupPara, setupTopo]

#######################################################################################################################
# Init Path
#######################################################################################################################
def initPath(nameFolder):
    basePath = pjoin(dirname(os.getcwd()), nameFolder)
    dataPath = pjoin(dirname(os.getcwd()), nameFolder, 'data')
    mdlPath  = pjoin(dirname(os.getcwd()), nameFolder, 'mdl')
    libPath  = pjoin(dirname(os.getcwd()), nameFolder, 'lib')
    resPath  = pjoin(dirname(os.getcwd()), nameFolder, 'results')
    parPath  = pjoin(dirname(os.getcwd()), nameFolder, 'para')
    setupPath = {'basePath': basePath, 'dataPath': dataPath, 'mdlPath': mdlPath, 'libPath': libPath, 'resPath': resPath, 'parPath': parPath}

    return setupPath

#######################################################################################################################
# Zero-Order-Hold Easy
#######################################################################################################################
def zoh_easy(x, c):
    xs = np.zeros(np.size(x))
    for i in range(0, int(len(x))):
        if i == 0:
            xs[i] = x[0]
            h_old = x[0]
        else:
            if c[i] > 0:
                xs[i] = x[i]
                h_old = x[i]
            elif c[i] < 0:
                xs[i] = 0
                h_old = 0
            else:
                xs[i] = h_old
    return xs
    
#######################################################################################################################
# Zero-Order Hold (ZOH)
#######################################################################################################################
def zoh(x, c, e, th):
    xs = np.zeros(np.size(x))
    hold = th
    if e == 'SE':
        for i in range(0, int(len(x))):
            if i == 0:
                xs[i] = x[0]
                i_old = i
            elif c[i] >= 0.99:
                if hold >=th:
                    xs[i] = x[i]
                    i_old = i
                    hold = 0
                else:
                    xs[i] = x[i_old]
                    hold = hold + 1
            else:
                xs[i] = x[i_old]
                hold = hold + 1
    else:
        for i in range(0, int(len(x))):
            if i == 0:
                xs[i] = x[0]
                i_old = i
            elif abs(c[i]) >= 0.99:
                if hold >=th:
                    xs[i] = x[i]
                    i_old = i
                    hold = 0
                else:
                    xs[i] = x[i_old]
                    hold = hold + 1
            else:
                xs[i] = x[i_old]
                hold = hold + 1

    return xs

#######################################################################################################################
# Convert to dB
#######################################################################################################################
def mag2dB(inp, ref):
    out = np.zeros(np.size(inp))
    for i in range(0, len(inp)):
        if inp[i] != 0:
            out[i] = 20 * np.log10(abs(inp[i]) / ref)
        else:
            out[i] = -60
            
    return out

#######################################################################################################################
# Order of Magnitude
#######################################################################################################################
def OoM(inp):
    return 10**math.floor(math.log(inp, 10))

#######################################################################################################################
# Root-Mean-Square
#######################################################################################################################
def rms(x):
    return np.sqrt(np.mean(x**2))

#######################################################################################################################
# Carrier Intersection
#######################################################################################################################
def cbInter(xs, c, Mi, tmin):
    if tmin == 0:
        if Mi != 0:
            s = xs > c
            s = s.astype(int) - (~s).astype(int)
        else:
            s = np.zeros(np.size(xs))
    else:
        hold = tmin
        s = np.zeros(np.size(xs))
        for i in range(0, len(xs)):
            if hold >= tmin:
                if Mi != 0:
                    temp = xs[i] > c[i]
                    s[i] = temp.astype(int) - (~temp).astype(int)
                hold = 0
            else:
                s[i] = s[i-1]
                hold = hold + 1
    return s

#######################################################################################################################
# sample
#######################################################################################################################
def con2dis(x, t, Ts):
    xs = np.zeros(np.size(x))
    k = 0
    for i in range(0, len(xs)):
        if t[i] >= k*Ts:
            xs[i] = x[i]
            k = k + 1
        else:
            xs[i] = xs[i-1]
    return xs

#######################################################################################################################
# THD
#######################################################################################################################
def thd(x, t, cf, K):
    dt = t[1] - t[0]
    T = t[-1] - t[0]
    N = len(x)
    X_eff = np.sqrt(1/T* np.sum(x**2 * dt))
    X_1 = cf*np.abs(fft(x)/N)[K]
    X_thd = np.sqrt(X_eff**2 - X_1**2)/X_1
    
    return X_thd

#######################################################################################################################
# WTHD
#######################################################################################################################
def wthd(x, t, cf, K):
    N = len(x)
    X_thd = 0
    X = cf*np.abs(fft(x)/N)
    for i in range(K+1, N):
        X_thd = X_thd + (X[i]/i)**2
    X_thd = np.sqrt(X_thd) 
    
    return X_thd

#######################################################################################################################
# Dead-time
#######################################################################################################################
def deadTime(s, Td):
    Nd = 0
    T1 = (s == 1)
    T2 = (s == -1)
    T1_out = np.zeros(np.size(T1))
    T2_out = np.zeros(np.size(T2))
    N = len(s)
    for i in range(2, N):
        if T1[i-1] == 0 and T1[i] == 1:
            Nd = Td
        if Nd > 0:
            T1_out[i] = 0
            Nd = Nd - 1
        else:
            T1_out[i] = T1[i]

        if T2[i-1] == 0 and T2[i] == 1:
            Nd = Td
        if Nd > 0:
            T2_out[i] = 0
            Nd = Nd - 1
        else:
            T2_out[i] = T2[i]
    
    s_out = T1_out - T2_out
    
    return s_out
