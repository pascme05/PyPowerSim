#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         plotStat_DAB
# Date:         05.02.2026
# Author:       Dr. Pascal A. Schirmer
# Version:      V.1.0
# Copyright:    Pascal Schirmer
#######################################################################################################################
#######################################################################################################################

#######################################################################################################################
# Import libs
#######################################################################################################################
# ==============================================================================
# Internal
# ==============================================================================
from src.general.helpFnc import OoM, thd

# ==============================================================================
# External
# ==============================================================================
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import matplotlib.gridspec as gridspec
from scipy.fft import fft
import matplotlib
matplotlib.use('TkAgg')


#######################################################################################################################
# Function
#######################################################################################################################
def plotStat_DAB(time, freq, setup):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("START: Plotting Stationary DAB Waveforms")

    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    timeSw = time['Sw']
    timeAc = time['Ac']
    timeDc = time['Dc']
    timeElec = time.get('Elec', {})
    timeLoss = time.get('Loss', {})
    timeTher = time.get('Ther', {})
    freqSw = freq['Sw']
    freqAc = freq['Ac']
    freqDc = freq['Dc']

    fel = setup['Top']['fel']
    fs = setup['Par']['PWM']['fs']
    fsim = setup['Exp']['fsim']
    Q = int(fs / fel)
    R = setup['Top']['R']
    L = setup['Top']['L']
    Mi = setup['Dat']['stat']['Mi']
    Vdc = setup['Dat']['stat']['Vdc']
    phiE = setup['Top']['phiE']
    n_tr = setup['Top'].get('n', 1)
    phiDAB = setup['Dat']['stat']['PhiDAB']
    down = int(setup['Dat']['stat']['cyc']) - 2
    down2 = int(fsim / fs / 100)
    if down2 < 1:
        down2 = 1

    t = timeSw['t'].values
    f = fsim * np.linspace(0, 0.5, int(len(t) / 2))

    ###################################################################################################################
    # Pre-Processing
    ###################################################################################################################
    K = int(np.round((t[-1] - t[0]) * fel))
    start = int((len(t) - 1) / K) * (K - 1)
    ende = len(t)

    # Slicing
    t = t[start:ende]
    if isinstance(timeSw, pd.DataFrame):
        timeSw = timeSw.iloc[start:ende]
    else:
        timeSw = timeSw[start:ende]

    for c1 in timeAc:
        if isinstance(timeAc[c1], pd.DataFrame):
            timeAc[c1] = timeAc[c1].iloc[start:ende]
        else:
            timeAc[c1] = timeAc[c1][start:ende]

    for c1 in timeDc:
        if isinstance(timeDc[c1], pd.DataFrame):
            timeDc[c1] = timeDc[c1].iloc[start:ende]
        else:
            timeDc[c1] = timeDc[c1][start:ende]

    # New: Slice loss and elec data too
    for c1 in timeElec:
        if isinstance(timeElec[c1], (pd.DataFrame, pd.Series)):
            timeElec[c1] = timeElec[c1].iloc[start:ende]
        elif isinstance(timeElec[c1], np.ndarray):
            timeElec[c1] = timeElec[c1][start:ende]
        elif isinstance(timeElec[c1], dict):
            for c2 in timeElec[c1]:
                if isinstance(timeElec[c1][c2], (pd.DataFrame, pd.Series)):
                    timeElec[c1][c2] = timeElec[c1][c2].iloc[start:ende]
                elif isinstance(timeElec[c1][c2], np.ndarray):
                    timeElec[c1][c2] = timeElec[c1][c2][start:ende]

    for c1 in timeLoss:
        if isinstance(timeLoss[c1], (pd.DataFrame, pd.Series)):
            timeLoss[c1] = timeLoss[c1].iloc[start:ende]
        elif isinstance(timeLoss[c1], np.ndarray):
            timeLoss[c1] = timeLoss[c1][start:ende]
        elif isinstance(timeLoss[c1], dict):
            for c2 in timeLoss[c1]:
                if isinstance(timeLoss[c1][c2], (pd.DataFrame, pd.Series)):
                    timeLoss[c1][c2] = timeLoss[c1][c2].iloc[start:ende]
                elif isinstance(timeLoss[c1][c2], np.ndarray):
                    timeLoss[c1][c2] = timeLoss[c1][c2][start:ende]

    for c1 in timeTher:
        if isinstance(timeTher[c1], (pd.DataFrame, pd.Series)):
            timeTher[c1] = timeTher[c1].iloc[start:ende]
        elif isinstance(timeTher[c1], np.ndarray):
            timeTher[c1] = timeTher[c1][start:ende]
        elif isinstance(timeTher[c1], dict):
            for c2 in timeTher[c1]:
                if isinstance(timeTher[c1][c2], (pd.DataFrame, pd.Series)):
                    timeTher[c1][c2] = timeTher[c1][c2].iloc[start:ende]
                elif isinstance(timeTher[c1][c2], np.ndarray):
                    timeTher[c1][c2] = timeTher[c1][c2][start:ende]

    # Time Axis Scaling
    t_plot = t - t[0]
    t_max = np.max(t_plot)
    if t_max < 1e-3:
        t_plot = t_plot * 1e6
        t_label = "time in ($\mu$s)"
    elif t_max < 1:
        t_plot = t_plot * 1e3
        t_label = "time in (ms)"
    else:
        t_label = "time in (s)"

    # Freq Axis Scaling
    f_max = fs
    f_max_plt = 200
    if f_max < 1e3:
        f_plot = f
        f_label = "$f$ in (Hz)"
    elif f_max < 1e6:
        f_plot = f / 1e3
        f_label = "$f$ in (kHz)"
    else:
        f_plot = f / 1e6
        f_label = "$f$ in (MHz)"

    # ==============================================================================
    # 1) Switching Function
    # ==============================================================================
    gs = gridspec.GridSpec(3, 2)
    pl.figure(figsize=(12, 10))
    txt = "Switching Functions DAB for PWM Control with: " \
          + "$V_{dc}$=" + str(Vdc) + "V, " \
          + "$M_{i}$=" + str(Mi) + "$ ,Q$=" + str(Q) \
          + ", Phi=" + f"{math.degrees(phiDAB):.2f}"
    pl.suptitle(txt, size=18)
    pl.subplots_adjust(hspace=0.35, wspace=0.20, left=0.075, right=0.925, top=0.90, bottom=0.075)

    # ------------------------------------------
    # Modulation
    # ------------------------------------------
    pl.subplot(gs[0, :])
    pl.plot(t_plot[::down2], timeSw['sA'][::down2])
    pl.plot(t_plot[::down2], timeSw['sC'][::down2])
    pl.title('Time-domain Switching Functions Primary vs. Secondary')
    pl.xlabel(t_label)
    pl.ylabel("$s_{x}(t)$ (p.u)")
    pl.legend(["$s_{a}$", "$s_{c}$"], loc='upper right')
    pl.grid('on')

    # ------------------------------------------
    # Switching Waveform
    # ------------------------------------------
    # Time-Domain
    pl.subplot(gs[1, 0])
    pl.plot(t_plot[::down2], timeSw['sA'][::down2])
    pl.plot(t_plot[::down2], timeSw['sB'][::down2])
    pl.title('Time-domain Switching Functions Primary')
    pl.xlabel(t_label)
    pl.ylabel("$s_{x}(t)$ (p.u)")
    pl.legend(["$s_{a}$", "$s_{b}$"], loc='upper right')
    pl.grid('on')

    # Freq-Domain
    pl.subplot(gs[2, 0])
    pl.stem(f_plot[::down][0:f_max_plt], freqSw['Sa'][::down][0:f_max_plt])
    pl.title('Frequency-domain Switching Function Primary')
    pl.xlabel(f_label)
    pl.ylabel("$S_{a}(f)$ (p.u)")
    pl.yscale('log')
    pl.ylim(0.0001, )
    pl.grid('on')

    # ------------------------------------------
    # Reference Waveform
    # ------------------------------------------
    # Time-Domain
    pl.subplot(gs[1, 1])
    pl.plot(t_plot[::down2], timeSw['sC'][::down2])
    pl.plot(t_plot[::down2], timeSw['sD'][::down2])
    pl.title('Time-domain Switching Functions Secondary')
    pl.xlabel(t_label)
    pl.ylabel("$s_{x}(t)$ (p.u)")
    pl.legend(["$s_{c}$", "$s_{d}$"], loc='upper right')
    pl.grid('on')

    # Freq-Domain
    pl.subplot(gs[2, 1])
    if 'Sc' in freqSw:
        pl.stem(f_plot[::down][0:f_max_plt], freqSw['Sc'][::down][0:f_max_plt])
    elif 'Sb' in freqSw:
        pl.stem(f_plot[::down][0:f_max_plt], freqSw['Sb'][::down][0:f_max_plt])
    pl.title('Frequency-domain Switching Function Secondary')
    pl.xlabel(f_label)
    pl.ylabel("$S_{c}(f)$ (p.u)")
    pl.yscale('log')
    pl.ylim(0.0001, )
    pl.grid('on')

    # ==============================================================================
    # 2) Primary Side Voltage / Current
    # ==============================================================================
    plt.figure()
    txt = "Currents and Voltages DAB Primary Bridge for PWM Control with: " \
          + "$V_{dc}$=" + str(Vdc) + "V, " \
          + "$M_{i}$=" + str(Mi) + "$ ,Q$=" + str(Q) \
          + ", Phi=" + f"{math.degrees(phiDAB):.2f}"
    plt.suptitle(txt, size=18)
    plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)

    # ------------------------------------------
    # Phase Voltage
    # ------------------------------------------
    # Time
    plt.subplot(4, 2, 1)
    plt.plot(t_plot[::down2], timeAc['v_ac_pri'][::down2])
    plt.plot(t_plot[::down2], timeAc['v_ac_sec_ref'][::down2])
    plt.plot(t_plot[::down2], timeAc['v_L'][::down2])
    plt.ylabel("$v_{ac,pri}(t)$ (V)")
    plt.title('Time-domain Voltage AC-Side Primary')
    plt.legend(["$v_{ac,pri}$", "$v_{ac,sec-ref}$", "$v_{ac,Lk}$"], loc='upper right')
    plt.xlabel(t_label)
    plt.grid('on')

    # Frequency
    plt.subplot(4, 2, 2)
    m1, s1, _ = plt.stem(f_plot[::down][0:f_max_plt], freqAc['v_ac_pri'][::down][0:f_max_plt])
    m2, s2, _ = plt.stem(f_plot[::down][0:f_max_plt], freqAc['v_L'][::down][0:f_max_plt])
    plt.setp(m1, color='C0')
    plt.setp(s1, color='C0')
    plt.setp(m2, color='C1')
    plt.setp(s2, color='C1')
    plt.ylabel("$V_{ac,pri}(f)$ (V)")
    plt.title('Frequency-domain Voltage AC-Side Primary')
    plt.xlabel(f_label)
    plt.yscale('log')
    plt.legend(handles=[m1, m2], labels=["$V_{ac,pri}$", "$V_{ac,Lk}$"], loc='upper right')
    plt.ylim(OoM(max(freqAc['v_ac_pri']))/1000, )
    plt.grid('on')

    # ------------------------------------------
    # Phase Current
    # ------------------------------------------
    # Time
    plt.subplot(4, 2, 3)
    plt.plot(t_plot[::down2], timeAc['i_ac_pri'][::down2])
    plt.ylabel("$i_{ac,pri}(t)$ (A)")
    plt.title('Time-domain Currents AC-Side Primary')
    plt.xlabel(t_label)
    plt.grid('on')

    # Frequency
    plt.subplot(4, 2, 4)
    plt.stem(f_plot[::down][0:f_max_plt], freqAc['i_ac_pri'][::down][0:f_max_plt])
    plt.ylabel("$I_{ac,pri}(f)$ (V)")
    plt.title('Frequency-domain Current AC-Side Primary')
    plt.xlabel(f_label)
    plt.yscale('log')
    plt.ylim(OoM(max(freqAc['i_ac_pri'])) / 1000, )
    plt.grid('on')

    # ------------------------------------------
    # DC-Link Voltage
    # ------------------------------------------
    # Time
    plt.subplot(4, 2, 5)
    plt.plot(t_plot[::down2], timeDc['v_dc_pri'][::down2], t_plot[::down2], np.mean(timeDc['v_dc_pri']) * np.ones(np.size(timeDc['i_dc_pri'][::down2])), '--')
    plt.ylabel("$v_{dc,pri}(t)$ (V)")
    plt.title('Time-domain Voltage DC-Side Primary')
    plt.xlabel(t_label)
    plt.legend(["$v_{dc,pri}$", "$V_{dc,pri,avg}$"], loc='upper right')
    plt.grid('on')

    # Frequency
    plt.subplot(4, 2, 6)
    plt.stem(f_plot[::down][0:f_max_plt], freqDc['v_dc_pri'][::down][0:f_max_plt])
    plt.ylabel("$V_{dc,pri}(f)$ (V)")
    plt.title('Frequency-domain Voltage DC-Side Primary')
    plt.xlabel(f_label)
    plt.yscale('log')
    plt.ylim(OoM(max(freqDc['v_dc_pri'])) / 1000000, )
    plt.grid('on')

    # ------------------------------------------
    # DC-Link Current
    # ------------------------------------------
    # Time
    plt.subplot(4, 2, 7)
    plt.plot(t_plot[::down2], timeDc['i_dc_pri'][::down2], color='tab:blue')
    plt.plot(t_plot[::down2], np.mean(timeDc['i_dc_pri']) * np.ones(np.size(timeDc['i_dc_pri'][::down2])), '--', color='tab:blue')
    plt.plot(t_plot[::down2], timeDc['i_c_pri'][::down2], color='tab:orange')
    plt.plot(t_plot[::down2], np.mean(timeDc['i_c_pri']) * np.ones(np.size(timeDc['i_dc_pri'][::down2])), '--', color='tab:orange')
    plt.plot(t_plot[::down2], timeDc['i_dc_in'][::down2], color='tab:green')
    plt.plot(t_plot[::down2], np.mean(timeDc['i_dc_in']) * np.ones(np.size(timeDc['i_dc_in'][::down2])), '--', color='tab:green')
    plt.ylabel("$i_{dc,pri}(t)$ (A)")
    plt.title('Time-domain Currents DC-Side Primary')
    plt.xlabel(t_label)
    plt.legend(["$i_{dc,pri}$", "$I_{dc,pri,avg}$", "$i_{c,pri}$", "$I_{c,pri,avg}$"], loc='upper right')
    plt.grid('on')

    # Frequency
    plt.subplot(4, 2, 8)
    m1, s1, _ = plt.stem(f_plot[::down][0:f_max_plt], freqDc['i_dc_pri'][::down][0:f_max_plt])
    m2, s2, _ = plt.stem(f_plot[::down][0:f_max_plt], freqDc['i_c_pri'][::down][0:f_max_plt])
    plt.setp(m1, color='C0')
    plt.setp(s1, color='C0')
    plt.setp(m2, color='C1')
    plt.setp(s2, color='C1')
    plt.ylabel("$I_{dc,pri}(f)$ (A)")
    plt.title('Frequency-domain Current DC-Side Primary')
    plt.xlabel(f_label)
    plt.yscale('log')
    plt.legend(handles=[m1, m2], labels=["$I_{dc,pri}$", "$I_{c,pri}$"], loc='upper right')
    plt.ylim(OoM(max(freqDc['i_c_pri'])) / 1000, )
    plt.grid('on')

    # ==============================================================================
    # 3) Secondary Side Voltage / Current
    # ==============================================================================
    plt.figure()
    txt = "Currents and Voltages DAB Secondary Bridge for PWM Control with: " \
          + "$V_{dc}$=" + str(Vdc) + "V, " \
          + "$M_{i}$=" + str(Mi) + "$ ,Q$=" + str(Q) \
          + ", Phi=" + f"{math.degrees(phiDAB):.2f}"
    plt.suptitle(txt, size=18)
    plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)

    # ------------------------------------------
    # Phase Voltage
    # ------------------------------------------
    # Time
    plt.subplot(4, 2, 1)
    plt.plot(t_plot[::down2], timeAc['v_ac_sec'][::down2])
    plt.ylabel("$v_{ac,sec}(t)$ (V)")
    plt.title('Time-domain Voltage AC-Side Secondary')
    plt.legend(["$v_{ac,sec}$"], loc='upper right')
    plt.xlabel(t_label)
    plt.grid('on')

    # Frequency
    plt.subplot(4, 2, 2)
    plt.stem(f_plot[::down][0:f_max_plt], freqAc['v_ac_sec'][::down][0:f_max_plt])
    plt.ylabel("$V_{ac,sec}(f)$ (V)")
    plt.title('Frequency-domain Voltage AC-Side Secondary')
    plt.xlabel(f_label)
    plt.yscale('log')
    plt.legend(labels=["$V_{ac,sec}$"], loc='upper right')
    plt.ylim(OoM(max(freqAc['v_ac_sec'])) / 1000, )
    plt.grid('on')

    # ------------------------------------------
    # Phase Current
    # ------------------------------------------
    # Time
    plt.subplot(4, 2, 3)
    plt.plot(t_plot[::down2], timeAc['i_ac_sec'][::down2])
    plt.ylabel("$i_{ac,sec}(t)$ (A)")
    plt.title('Time-domain Currents AC-Side Secondary')
    plt.xlabel(t_label)
    plt.grid('on')

    # Frequency
    plt.subplot(4, 2, 4)
    plt.stem(f_plot[::down][0:f_max_plt], freqAc['i_ac_sec'][::down][0:f_max_plt])
    plt.ylabel("$I_{ac,sec}(f)$ (A)")
    plt.title('Frequency-domain Current AC-Side Secondary')
    plt.xlabel(f_label)
    plt.yscale('log')
    plt.ylim(OoM(max(freqAc['i_ac_sec'])) / 1000, )
    plt.grid('on')

    # ------------------------------------------
    # DC-Link Voltage
    # ------------------------------------------
    # Time
    plt.subplot(4, 2, 5)
    plt.plot(t_plot[::down2], timeDc['v_dc_sec'][::down2], t_plot[::down2],
             np.mean(timeDc['v_dc_sec']) * np.ones(np.size(timeDc['i_dc_sec'][::down2])), '--')
    plt.ylabel("$v_{dc,sec}(t)$ (V)")
    plt.title('Time-domain Voltage DC-Side Secondary')
    plt.xlabel(t_label)
    plt.legend(["$v_{dc,sec}$", "$V_{dc,sec,avg}$"], loc='upper right')
    plt.grid('on')

    # Frequency
    plt.subplot(4, 2, 6)
    plt.stem(f_plot[::down][0:f_max_plt], freqDc['v_dc_sec'][::down][0:f_max_plt])
    plt.ylabel("$V_{dc,sec}(f)$ (V)")
    plt.title('Frequency-domain Voltage DC-Side Secondary')
    plt.xlabel(f_label)
    plt.yscale('log')
    plt.ylim(OoM(max(freqDc['v_dc_sec'])) / 1000000, )
    plt.grid('on')

    # ------------------------------------------
    # DC-Link Current
    # ------------------------------------------
    # Time
    plt.subplot(4, 2, 7)
    plt.plot(t_plot[::down2], timeDc['i_dc_sec'][::down2], color='tab:blue')
    plt.plot(t_plot[::down2], np.mean(timeDc['i_dc_sec']) * np.ones(np.size(timeDc['i_dc_sec'][::down2])), '--', color='tab:blue')
    plt.plot(t_plot[::down2], timeDc['i_c_sec'][::down2], color='tab:orange')
    plt.plot(t_plot[::down2], np.mean(timeDc['i_c_sec']) * np.ones(np.size(timeDc['i_c_sec'][::down2])), '--', color='tab:orange')
    plt.plot(t_plot[::down2], timeDc['i_dc_out'][::down2], color='tab:green')
    plt.plot(t_plot[::down2], np.mean(timeDc['i_dc_out']) * np.ones(np.size(timeDc['i_dc_out'][::down2])), '--', color='tab:green')
    plt.ylabel("$i_{dc,sec}(t)$ (A)")
    plt.title('Time-domain Currents DC-Side Secondary')
    plt.xlabel(t_label)
    plt.legend(["$i_{dc,sec}$", "$I_{dc,pri,sec}$", "$i_{c,sec}$", "$I_{c,pri,sec}$"], loc='upper right')
    plt.grid('on')

    # Frequency
    plt.subplot(4, 2, 8)
    m1, s1, _ = plt.stem(f_plot[::down][0:f_max_plt], freqDc['i_dc_sec'][::down][0:f_max_plt])
    m2, s2, _ = plt.stem(f_plot[::down][0:f_max_plt], freqDc['i_c_sec'][::down][0:f_max_plt])
    plt.setp(m1, color='C0')
    plt.setp(s1, color='C0')
    plt.setp(m2, color='C1')
    plt.setp(s2, color='C1')
    plt.ylabel("$I_{dc,sec}(f)$ (A)")
    plt.title('Frequency-domain Current DC-Side Secondary')
    plt.xlabel(f_label)
    plt.yscale('log')
    plt.legend(handles=[m1, m2], labels=["$I_{dc,sec}$", "$I_{c,sec}$"], loc='upper right')
    plt.ylim(OoM(max(freqDc['i_c_sec'])) / 1000, )
    plt.grid('on')

    # ==============================================================================
    # 4) Secondary Side Current
    # ==============================================================================
    plt.figure()
    txt = "Currents DAB Secondary Bridge for PWM Control with: " + "$V_{dc}$=" + str(Vdc) + "V, " + "$M_{i}$=" + str(
        Mi) + "$ ,Q$=" + str(Q) + ", Phi=" + str(phiDAB)
    plt.suptitle(txt, size=18)
    plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)

    # ------------------------------------------
    # Phase
    # ------------------------------------------
    # Time
    plt.subplot(2, 2, 1)
    plt.plot(t_plot[::down2], timeAc['i_ac_sec'][::down2])
    plt.ylabel("$i_{s}(t)$ (A)")
    plt.title('Time-domain Currents AC-Side Secondary')
    plt.xlabel(t_label)
    plt.grid('on')

    # Frequency
    ax = plt.subplot(2, 2, 2)
    i_s_spec = freqAc.get('i_ac_sec', freqAc.get('I_b'))
    if i_s_spec is not None:
        plt.stem(f_plot[::down][0:50], i_s_spec[::down][0:50])
    plt.ylabel("$I_{s}(f)$ (A)")
    plt.xlim(0, 50)
    plt.title('Frequency-domain Current AC-Side Secondary')
    plt.xlabel("$f/f_{1}$ (Hz/Hz)")
    plt.yscale('log')
    if i_s_spec is not None:
        plt.ylim(0.1 / OoM(max(i_s_spec)), )
    plt.grid('on')

    # ------------------------------------------
    # DC-Link
    # ------------------------------------------
    # Time
    plt.subplot(2, 2, 3)
    plt.plot(t_plot[::down2], timeDc['i_dc_sec'][::down2], t_plot[::down2], np.mean(timeDc['i_dc_sec']) * np.ones(np.size(timeDc['i_dc_sec'][::down2])), '--')
    plt.ylabel("$i_{dc,sec}(t)$ (A)")
    plt.title('Time-domain Currents DC-Side Secondary')
    plt.xlabel(t_label)
    plt.legend(["$i_{dc,sec}$", "$I_{dc,sec,avg}$"], loc='upper right')
    plt.grid('on')

    # Frequency
    ax = plt.subplot(2, 2, 4)
    i_dc_s_spec = freqDc.get('i_dc_sec', freqDc.get('I_dc_s'))
    if i_dc_s_spec is not None:
        plt.stem(f_plot[::down][0:50], i_dc_s_spec[::down][0:50])
    plt.ylabel("$I_{dc,sec}(f)$ (A)")
    plt.xlim(0, 50)
    plt.title('Frequency-domain Current DC-Side Secondary')
    plt.xlabel("$f/f_{1}$ (Hz/Hz)")
    plt.yscale('log')
    if i_dc_s_spec is not None:
        plt.ylim(0.1 / OoM(max(i_dc_s_spec)), )
    plt.grid('on')

    # ==============================================================================
    # 5) Secondary Side Voltage
    # ==============================================================================
    plt.figure()
    txt = "Voltages DAB Secondary Bridge for PWM Control with: " + "$V_{dc}$=" + str(Vdc) + "V, " + "$M_{i}$=" + str(
        Mi) + "$ ,Q$=" + str(Q) + ", Phi=" + str(phiDAB)
    plt.suptitle(txt, size=18)
    plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)

    # ------------------------------------------
    # Phase
    # ------------------------------------------
    # Time
    plt.subplot(2, 2, 1)
    plt.plot(t_plot[::down2], timeAc['v_ac_sec'][::down2])
    plt.ylabel("$v_{s}(t)$ (V)")
    plt.title('Time-domain Voltages AC-Side Secondary')
    plt.xlabel(t_label)
    plt.grid('on')

    # Frequency
    ax = plt.subplot(2, 2, 2)
    v_s_spec = freqAc.get('v_ac_sec', freqAc.get('V_b'))
    if v_s_spec is None and 'v_b0' in freqAc: v_s_spec = freqAc['v_b0']
    if v_s_spec is not None:
        plt.stem(f_plot[::down][0:50], v_s_spec[::down][0:50])
    plt.ylabel("$V_{s}(f)$ (V)")
    plt.xlim(0, 50)
    plt.title('Frequency-domain Voltages AC-Side Secondary')
    plt.xlabel("$f/f_{1}$ (Hz/Hz)")
    plt.yscale('log')
    if v_s_spec is not None:
        plt.ylim(0.1 / OoM(max(v_s_spec)), )
    plt.grid('on')

    # ------------------------------------------
    # DC-Link
    # ------------------------------------------
    # Time
    plt.subplot(2, 2, 3)
    plt.plot(t_plot[::down2], timeDc['v_dc_sec'][::down2], t_plot[::down2], np.mean(timeDc['v_dc_sec']) * np.ones(np.size(timeDc['v_dc_sec'][::down2])), '--')
    plt.ylabel("$v_{dc,sec}(t)$ (V)")
    plt.title('Time-domain Voltages DC-Side Secondary')
    plt.xlabel(t_label)
    plt.legend(["$v_{dc,sec}$", "$V_{dc,sec,avg}$"], loc='upper right')
    plt.grid('on')

    # Frequency
    ax = plt.subplot(2, 2, 4)
    v_dc_s_spec = freqDc.get('v_dc_sec', freqDc.get('V_dc_s'))
    if v_dc_s_spec is not None:
        plt.stem(f_plot[::down][0:50], v_dc_s_spec[::down][0:50])
    plt.ylabel("$V_{dc,sec}(f)$ (V)")
    plt.xlim(0, 50)
    plt.title('Frequency-domain Voltages DC-Side Secondary')
    plt.xlabel("$f/f_{1}$ (Hz/Hz)")
    plt.yscale('log')
    if v_dc_s_spec is not None:
        plt.ylim(0.1 / OoM(max(v_dc_s_spec)), )
    plt.grid('on')

    # ==============================================================================
    # 6) Time-domain Switches
    # ==============================================================================
    if 'sw' in timeElec:
        # IDs
        id_sw = ['S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8']
        id_T = ['T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8']
        id_D = ['D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8']
        id_C = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8']
        Ta = setup['Dat']['stat']['Tc']

        plt.figure(figsize=(18, 12))
        txt = "Time domain switches DAB for PWM control with: " + "$V_{dc}$=" + str(Vdc) + "V, " + "$M_{i}$=" + str(
            Mi) + "$ ,Q$=" + str(Q) + ", Phi=" + str(phiDAB)
        plt.suptitle(txt, size=18)
        plt.subplots_adjust(hspace=0.45, wspace=0.35, left=0.05, right=0.95, top=0.90, bottom=0.075)

        for i in range(0, 8):
            sid = id_sw[i]
            if sid in timeElec['sw']:
                # Current
                plt.subplot(8, 5, 5 * i + 1)
                plt.plot(t_plot[::down2], timeElec['sw'][sid]['i_T'][::down2], t_plot[::down2], timeElec['sw'][sid]['i_D'][::down2])
                plt.ylabel("$i(t)$ (A)")
                plt.title('Currents ' + sid)
                if i == 7: plt.xlabel(t_label)
                else: plt.xticks([], [])
                plt.grid('on')

                # Voltage
                plt.subplot(8, 5, 5 * i + 2)
                plt.plot(t_plot[::down2], timeElec['sw'][sid]['v_T'][::down2], t_plot[::down2], timeElec['sw'][sid]['v_D'][::down2])
                plt.ylabel("$v(t)$ (V)")
                plt.title('Voltages ' + sid)
                if i == 7: plt.xlabel(t_label)
                else: plt.xticks([], [])
                plt.grid('on')

                # Conduction Losses
                plt.subplot(8, 5, 5 * i + 3)
                if sid in timeLoss['sw']:
                    plt.plot(t_plot[::down2], timeLoss['sw'][sid]['p_T_c'][::down2], t_plot[::down2], timeLoss['sw'][sid]['p_D_c'][::down2])
                plt.ylabel("$p_{c}(t)$ (W)")
                plt.title('Cond. Losses ' + sid)
                if i == 7: plt.xlabel(t_label)
                else: plt.xticks([], [])
                plt.grid('on')

                # Switching Losses
                plt.subplot(8, 5, 5 * i + 4)
                if sid in timeLoss['sw']:
                    plt.plot(t_plot[::down2], timeLoss['sw'][sid]['p_T_s'][::down2], t_plot[::down2], timeLoss['sw'][sid]['p_D_s'][::down2])
                plt.ylabel("$p_{s}(t)$ (W)")
                plt.title('Swi. Losses ' + sid)
                if i == 7: plt.xlabel(t_label)
                else: plt.xticks([], [])
                plt.grid('on')

                # Temperature
                plt.subplot(8, 5, 5 * i + 5)
                if 'sw' in timeTher:
                    plt.plot(t_plot[::down2], timeTher['sw'][id_T[i]][::down2], t_plot[::down2], timeTher['sw'][id_D[i]][::down2], 
                             t_plot[::down2], timeTher['sw'][id_C[i]][::down2], t_plot[::down2], Ta * np.ones(np.size(t_plot[::down2])))
                plt.ylabel("$\Theta(t)$ (°C)")
                plt.title('Thermal ' + sid)
                if i == 7: plt.xlabel(t_label)
                else: plt.xticks([], [])
                plt.grid('on')

    # ==============================================================================
    # 7) Time-domain Capacitor
    # ==============================================================================
    if 'cap' in timeElec:
        plt.figure(figsize=(12, 10))
        plt.suptitle("Time domain capacitor DAB", size=18)
        plt.subplots_adjust(hspace=0.45, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)

        caps = list(timeElec['cap'].keys())
        for i, cid in enumerate(caps):
            # Current
            plt.subplot(len(caps), 4, 4 * i + 1)
            plt.plot(t_plot[::down2], timeElec['cap'][cid]['i_c'][::down2])
            plt.ylabel("$i(t)$ (A)")
            plt.title('Current ' + cid)
            plt.grid('on')

            # Voltage
            plt.subplot(len(caps), 4, 4 * i + 2)
            plt.plot(t_plot[::down2], timeElec['cap'][cid]['v_c'][::down2])
            plt.ylabel("$v(t)$ (V)")
            plt.title('Voltage ' + cid)
            plt.grid('on')

            # Losses
            plt.subplot(len(caps), 4, 4 * i + 3)
            if 'cap' in timeLoss and cid in timeLoss['cap']:
                plt.plot(t_plot[::down2], timeLoss['cap'][cid]['p_L'][::down2])
            plt.ylabel("$p(t)$ (W)")
            plt.title('Losses ' + cid)
            plt.grid('on')

            # Temperature
            plt.subplot(len(caps), 4, 4 * i + 4)
            if 'cap' in timeTher and cid in timeTher['cap']:
                plt.plot(t_plot[::down2], timeTher['cap'][cid][::down2])
            plt.ylabel("$\Theta(t)$ (°C)")
            plt.title('Thermal ' + cid)
            plt.grid('on')

    # ==============================================================================
    # 8) Transformer Overview
    # ==============================================================================
    if 'tra' in timeElec and len(timeElec['tra']) > 0:
        plt.figure(figsize=(12, 10))
        plt.suptitle("Transformer Overview: Electrical and Thermal", size=18)
        plt.subplots_adjust(hspace=0.5, wspace=0.3)

        # Electrical: Flux and B
        plt.subplot(3, 2, 1)
        plt.plot(t_plot[::down2], timeElec['tra']['Phi'][::down2])
        plt.title('Magnetic Flux')
        plt.xlabel(t_label)
        plt.ylabel('$\Phi$ (Vs)')
        plt.grid('on')

        plt.subplot(3, 2, 2)
        plt.plot(t_plot[::down2], timeElec['tra']['B'][::down2])
        plt.title('Magnetic Flux Density')
        plt.xlabel(t_label)
        plt.ylabel('B (T)')
        plt.grid('on')

        # Electrical: Magnetizing Current
        plt.subplot(3, 2, 3)
        plt.plot(t_plot[::down2], timeElec['tra']['i_m'][::down2])
        plt.title('Magnetizing Current')
        plt.xlabel(t_label)
        plt.ylabel('$i_{m}$ (A)')
        plt.grid('on')

        # Losses: Overview
        plt.subplot(3, 2, 4)
        if 'tra' in timeLoss:
            p_pc = np.mean(timeLoss['tra']['p_PC'])
            p_sc = np.mean(timeLoss['tra']['p_SC'])
            p_cc = np.mean(timeLoss['tra']['p_CC'])
            plt.bar(['Pri Cu', 'Sec Cu', 'Core'], [p_pc, p_sc, p_cc])
            plt.title('Average Transformer Losses')
            plt.ylabel('Losses (W)')
            plt.grid('on', axis='y')

        # Thermal Overview
        plt.subplot(3, 1, 3)
        if 'tra' in timeTher:
            for nid in timeTher['tra']:
                plt.plot(t_plot[::down2], timeTher['tra'][nid][::down2], label=nid)
        plt.title('Transformer Temperatures')
        plt.xlabel(t_label)
        plt.ylabel('Temperature (°C)')
        plt.legend()
        plt.grid('on')

    plt.show()

    ###################################################################################################################
    # MSG OUT
    ###################################################################################################################
    print("END: Plotting Stationary DAB Waveforms")
