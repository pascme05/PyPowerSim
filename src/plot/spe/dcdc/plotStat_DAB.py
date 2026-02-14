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
    f_scale = 10000
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
    ax1 = plt.subplot(4, 2, 1)
    plt.plot(t_plot[::down2], timeAc['v_ac_pri'][::down2])
    plt.plot(t_plot[::down2], timeAc['v_L'][::down2])
    plt.plot(t_plot[::down2], timeAc['v_ac_sec_ref'][::down2])
    plt.ylabel("$v_{ac,pri}(t)$ (V)")
    plt.title('Time-domain Voltage AC-Side Primary')
    plt.legend(["$v_{ac,pri}$", "$v_{ac,Lk}$", "$v_{ac,sec-ref}$"], loc='upper right')
    # plt.xlabel(t_label)
    plt.grid('on')

    # Frequency
    ax2 = plt.subplot(4, 2, 2)
    m1, s1, _ = plt.stem(f_plot[::down][0:f_max_plt], freqAc['v_ac_pri'][::down][0:f_max_plt])
    m2, s2, _ = plt.stem(f_plot[::down][0:f_max_plt], freqAc['v_L'][::down][0:f_max_plt])
    plt.setp(m1, color='C0')
    plt.setp(s1, color='C0')
    plt.setp(m2, color='C1')
    plt.setp(s2, color='C1')
    plt.ylabel("$V_{ac,pri}(f)$ (V)")
    plt.title('Frequency-domain Voltage AC-Side Primary')
    # plt.xlabel(f_label)
    plt.yscale('log')
    plt.legend(handles=[m1, m2], labels=["$V_{ac,pri}$", "$V_{ac,Lk}$"], loc='upper right')
    plt.ylim(OoM(max(freqAc['v_ac_pri']))/f_scale, )
    plt.grid('on')

    # ------------------------------------------
    # Phase Current
    # ------------------------------------------
    # Time
    plt.subplot(4, 2, 3, sharex=ax1)
    plt.plot(t_plot[::down2], timeAc['i_ac_pri'][::down2])
    plt.ylabel("$i_{ac,pri}(t)$ (A)")
    plt.title('Time-domain Currents AC-Side Primary')
    # plt.xlabel(t_label)
    plt.grid('on')

    # Frequency
    plt.subplot(4, 2, 4, sharex=ax2)
    plt.stem(f_plot[::down][0:f_max_plt], freqAc['i_ac_pri'][::down][0:f_max_plt])
    plt.ylabel("$I_{ac,pri}(f)$ (A)")
    plt.title('Frequency-domain Current AC-Side Primary')
    # plt.xlabel(f_label)
    plt.yscale('log')
    plt.ylim(OoM(max(freqAc['i_ac_pri'])) / f_scale, )
    plt.grid('on')

    # ------------------------------------------
    # DC-Link Voltage
    # ------------------------------------------
    # Time
    plt.subplot(4, 2, 5, sharex=ax1)
    plt.plot(t_plot[::down2], timeDc['v_dc_pri'][::down2], t_plot[::down2], np.mean(timeDc['v_dc_pri']) * np.ones(np.size(timeDc['i_dc_pri'][::down2])))
    plt.ylabel("$v_{dc,pri}(t)$ (V)")
    plt.title('Time-domain Voltage DC-Side Primary')
    # plt.xlabel(t_label)
    plt.legend(["$v_{dc,pri}$", "$V_{dc,pri,avg}$"], loc='upper right')
    plt.grid('on')

    # Frequency
    plt.subplot(4, 2, 6, sharex=ax2)
    plt.stem(f_plot[::down][0:f_max_plt], freqDc['v_dc_pri'][::down][0:f_max_plt])
    plt.ylabel("$V_{dc,pri}(f)$ (V)")
    plt.title('Frequency-domain Voltage DC-Side Primary')
    # plt.xlabel(f_label)
    plt.yscale('log')
    plt.ylim(OoM(max(freqDc['v_dc_pri'])) / 1000000, )
    plt.grid('on')

    # ------------------------------------------
    # DC-Link Current
    # ------------------------------------------
    # Time
    plt.subplot(4, 2, 7, sharex=ax1)
    plt.plot(t_plot[::down2], timeDc['i_dc_pri'][::down2], color='tab:blue')
    plt.plot(t_plot[::down2], timeDc['i_c_pri'][::down2], color='tab:orange')
    plt.plot(t_plot[::down2], timeDc['i_dc_in'][::down2], color='tab:green')
    plt.ylabel("$i_{dc,pri}(t)$ (A)")
    plt.title('Time-domain Currents DC-Side Primary')
    plt.xlabel(t_label)
    plt.legend(["$i_{dc,pri}$", "$i_{c,pri}$", "$I_{dc,in}$"], loc='upper right')
    plt.grid('on')

    # Frequency
    plt.subplot(4, 2, 8, sharex=ax2)
    m1, s1, _ = plt.stem(f_plot[::down][0:f_max_plt], freqDc['i_dc_pri'][::down][0:f_max_plt])
    m2, s3, _ = plt.stem(f_plot[::down][0:f_max_plt], freqDc['i_c_pri'][::down][0:f_max_plt])
    m3, s3, _ = plt.stem(f_plot[::down][0:f_max_plt], freqDc['i_dc_in'][::down][0:f_max_plt])
    plt.setp(m1, color='C0')
    plt.setp(s1, color='C0')
    plt.setp(m2, color='C1')
    plt.setp(s2, color='C1')
    plt.setp(m3, color='C2')
    plt.setp(s3, color='C2')
    plt.ylabel("$I_{dc,pri}(f)$ (A)")
    plt.title('Frequency-domain Current DC-Side Primary')
    plt.xlabel(f_label)
    plt.yscale('log')
    plt.legend(handles=[m1, m2, m3], labels=["$I_{dc,pri}$", "$I_{c,pri}$", "$I_{dc,in}$"], loc='upper right')
    plt.ylim(OoM(max(freqDc['i_dc_pri'])) / f_scale, )
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
    ax3 = plt.subplot(4, 2, 1)
    plt.plot(t_plot[::down2], timeAc['v_ac_sec'][::down2])
    plt.ylabel("$v_{ac,sec}(t)$ (V)")
    plt.title('Time-domain Voltage AC-Side Secondary')
    plt.legend(["$v_{ac,sec}$"], loc='upper right')
    # plt.xlabel(t_label)
    plt.grid('on')

    # Frequency
    ax4 = plt.subplot(4, 2, 2)
    plt.stem(f_plot[::down][0:f_max_plt], freqAc['v_ac_sec'][::down][0:f_max_plt])
    plt.ylabel("$V_{ac,sec}(f)$ (V)")
    plt.title('Frequency-domain Voltage AC-Side Secondary')
    # plt.xlabel(f_label)
    plt.yscale('log')
    plt.legend(labels=["$V_{ac,sec}$"], loc='upper right')
    plt.ylim(OoM(max(freqAc['v_ac_sec'])) / f_scale, )
    plt.grid('on')

    # ------------------------------------------
    # Phase Current
    # ------------------------------------------
    # Time
    plt.subplot(4, 2, 3, sharex=ax3)
    plt.plot(t_plot[::down2], timeAc['i_ac_sec'][::down2])
    plt.ylabel("$i_{ac,sec}(t)$ (A)")
    plt.title('Time-domain Currents AC-Side Secondary')
    # plt.xlabel(t_label)
    plt.grid('on')

    # Frequency
    plt.subplot(4, 2, 4, sharex=ax4)
    plt.stem(f_plot[::down][0:f_max_plt], freqAc['i_ac_sec'][::down][0:f_max_plt])
    plt.ylabel("$I_{ac,sec}(f)$ (A)")
    plt.title('Frequency-domain Current AC-Side Secondary')
    # plt.xlabel(f_label)
    plt.yscale('log')
    plt.ylim(OoM(max(freqAc['i_ac_sec'])) / f_scale, )
    plt.grid('on')

    # ------------------------------------------
    # DC-Link Voltage
    # ------------------------------------------
    # Time
    plt.subplot(4, 2, 5, sharex=ax3)
    plt.plot(t_plot[::down2], timeDc['v_dc_sec'][::down2], t_plot[::down2],
             np.mean(timeDc['v_dc_sec']) * np.ones(np.size(timeDc['i_dc_sec'][::down2])))
    plt.ylabel("$v_{dc,sec}(t)$ (V)")
    plt.title('Time-domain Voltage DC-Side Secondary')
    # plt.xlabel(t_label)
    plt.legend(["$v_{dc,sec}$", "$V_{dc,sec,avg}$"], loc='upper right')
    plt.grid('on')

    # Frequency
    plt.subplot(4, 2, 6, sharex=ax4)
    plt.stem(f_plot[::down][0:f_max_plt], freqDc['v_dc_sec'][::down][0:f_max_plt])
    plt.ylabel("$V_{dc,sec}(f)$ (V)")
    plt.title('Frequency-domain Voltage DC-Side Secondary')
    # plt.xlabel(f_label)
    plt.yscale('log')
    plt.ylim(OoM(max(freqDc['v_dc_sec'])) / 1000000, )
    plt.grid('on')

    # ------------------------------------------
    # DC-Link Current
    # ------------------------------------------
    # Time
    plt.subplot(4, 2, 7, sharex=ax3)
    plt.plot(t_plot[::down2], timeDc['i_dc_sec'][::down2], color='tab:blue')
    plt.plot(t_plot[::down2], timeDc['i_c_sec'][::down2], color='tab:orange')
    plt.plot(t_plot[::down2], timeDc['i_dc_out'][::down2], color='tab:green')
    plt.ylabel("$i_{dc,sec}(t)$ (A)")
    plt.title('Time-domain Currents DC-Side Secondary')
    plt.xlabel(t_label)
    plt.legend(["$i_{dc,sec}$", "$i_{c,sec}$", "$i_{dc,out}$"], loc='upper right')
    plt.grid('on')

    # Frequency
    plt.subplot(4, 2, 8, sharex=ax4)
    m1, s1, _ = plt.stem(f_plot[::down][0:f_max_plt], freqDc['i_dc_sec'][::down][0:f_max_plt])
    m2, s2, _ = plt.stem(f_plot[::down][0:f_max_plt], freqDc['i_c_sec'][::down][0:f_max_plt])
    m3, s3, _ = plt.stem(f_plot[::down][0:f_max_plt], freqDc['i_dc_out'][::down][0:f_max_plt])
    plt.setp(m1, color='C0')
    plt.setp(s1, color='C0')
    plt.setp(m2, color='C1')
    plt.setp(s2, color='C1')
    plt.setp(m3, color='C2')
    plt.setp(s3, color='C2')
    plt.ylabel("$I_{dc,sec}(f)$ (A)")
    plt.title('Frequency-domain Current DC-Side Secondary')
    plt.xlabel(f_label)
    plt.yscale('log')
    plt.legend(handles=[m1, m2, m3], labels=["$I_{dc,sec}$", "$I_{c,sec}$", "$I_{dc,out}$"], loc='upper right')
    plt.ylim(OoM(max(freqDc['i_dc_sec'])) / f_scale, )
    plt.grid('on')

    # ==============================================================================
    # 4) Time-domain Switches
    # ==============================================================================
    if 'sw' in timeElec:
        # IDs
        id_sw = ['S1', 'S2', 'S5', 'S6']
        id_T = ['T1', 'T2', 'T5', 'T6']
        id_D = ['D1', 'D2', 'D5', 'D6']
        id_C = ['C1', 'C2', 'C5', 'C6']
        Ta = setup['Dat']['stat']['Tc']

        plt.figure(figsize=(18, 12))
        txt = "Time domain switches DAB for PWM control with: " \
              + "$V_{dc}$=" + str(Vdc) + "V, " \
              + "$M_{i}$=" + str(Mi) + "$ ,Q$=" + str(Q) \
              + ", Phi=" + f"{math.degrees(phiDAB):.2f}"
        plt.suptitle(txt, size=18)
        plt.subplots_adjust(hspace=0.45, wspace=0.35, left=0.05, right=0.95, top=0.90, bottom=0.075)

        ax_list = []
        for i in range(0, len(id_sw)):
            sid = id_sw[i]
            if sid in timeElec['sw']:
                # Current
                if i == 0:
                    ax_first = plt.subplot(5, len(id_sw), i+1)
                else:
                    ax_first = plt.subplot(5, len(id_sw), i+1, sharex=ax_first)
                
                plt.plot(t_plot[::down2], timeElec['sw'][sid]['i_T'][::down2], t_plot[::down2], timeElec['sw'][sid]['i_D'][::down2])
                plt.ylabel("$i(t)$ (A)")
                plt.title('Currents ' + sid)
                plt.tick_params(axis='x', which='both', length=0, labelbottom=False)
                plt.grid('on')
                plt.legend(["$Trans$", "$Diode$"], loc='upper right')

                # Voltage
                plt.subplot(5, len(id_sw), i+5, sharex=ax_first)
                plt.plot(t_plot[::down2], timeElec['sw'][sid]['v_T'][::down2], t_plot[::down2], timeElec['sw'][sid]['v_D'][::down2])
                plt.ylabel("$v(t)$ (V)")
                plt.title('Voltages ' + sid)
                # plt.xticks([], [])
                plt.tick_params(axis='x', which='both', length=0, labelbottom=False)
                plt.grid('on')
                plt.legend(["$Trans$", "$Diode$"], loc='upper right')

                # Conduction Losses
                plt.subplot(5, len(id_sw), i+9, sharex=ax_first)
                if sid in timeLoss['sw']:
                    plt.plot(t_plot[::down2], timeLoss['sw'][sid]['p_T_c'][::down2], t_plot[::down2], timeLoss['sw'][sid]['p_D_c'][::down2])
                plt.ylabel("$p_{c}(t)$ (W)")
                plt.title('Conduction Losses ' + sid)
                plt.tick_params(axis='x', which='both', length=0, labelbottom=False)
                plt.grid('on')
                plt.legend(["$Trans$", "$Diode$"], loc='upper right')

                # Switching Losses
                plt.subplot(5, len(id_sw), i+13, sharex=ax_first)
                if sid in timeLoss['sw']:
                    plt.plot(t_plot[::down2], timeLoss['sw'][sid]['p_T_s'][::down2], t_plot[::down2], timeLoss['sw'][sid]['p_D_s'][::down2])
                plt.ylabel("$p_{s}(t)$ (W)")
                plt.title('Switching Losses ' + sid)
                plt.tick_params(axis='x', which='both', length=0, labelbottom=False)
                plt.grid('on')
                plt.legend(["$Trans$", "$Diode$"], loc='upper right')

                # Temperature
                plt.subplot(5, len(id_sw), i+17, sharex=ax_first)
                if 'sw' in timeTher:
                    plt.plot(t_plot[::down2], timeTher['sw'][id_T[i]][::down2], t_plot[::down2], timeTher['sw'][id_D[i]][::down2], 
                             t_plot[::down2], timeTher['sw'][id_C[i]][::down2], t_plot[::down2], Ta * np.ones(np.size(t_plot[::down2])))
                plt.ylabel("$\Theta(t)$ (°C)")
                plt.title('Thermal ' + sid)
                plt.xlabel(t_label)
                plt.grid('on')
                plt.legend(["$Trans$", "$Diode$", "Case", "Ambient"], loc='upper right')

    # ==============================================================================
    # 5) Time-domain Capacitor
    # ==============================================================================
    if 'cap' in timeElec:
        plt.figure(figsize=(12, 10))
        txt = "Time domain capacitor DAB for PWM control with: " \
              + "$V_{dc}$=" + str(Vdc) + "V, " \
              + "$M_{i}$=" + str(Mi) + "$ ,Q$=" + str(Q) \
              + ", Phi=" + f"{math.degrees(phiDAB):.2f}"
        plt.suptitle(txt, size=18)
        plt.subplots_adjust(hspace=0.45, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)

        caps = list(timeElec['cap'].keys())
        n_caps = max(len(caps), 1)

        # Two-column layout (one column per capacitor):
        # rows: current, voltage, losses, temperature
        for i, cid in enumerate(caps):
            # Current (row 0)
            ax_i = plt.subplot(4, n_caps, 0 * n_caps + i + 1)
            plt.plot(t_plot[::down2], timeElec['cap'][cid]['i_c'][::down2])
            plt.plot(t_plot[::down2], np.mean(timeElec['cap'][cid]['i_c'][::down2]) * np.ones(np.size(t_plot[::down2])))
            plt.ylabel("$i(t)$ (A)")
            plt.title(f"Current {cid}")
            plt.legend(["$i_{c}$", "$I_{c,avg}$"], loc='upper right')
            plt.grid('on')

            # Voltage (row 1)
            ax_v = plt.subplot(4, n_caps, 1 * n_caps + i + 1, sharex=ax_i)
            plt.plot(t_plot[::down2], timeElec['cap'][cid]['v_c'][::down2])
            plt.plot(t_plot[::down2], np.mean(timeElec['cap'][cid]['v_c'][::down2]) * np.ones(np.size(t_plot[::down2])))
            plt.ylabel("$v(t)$ (V)")
            plt.title(f"Voltage {cid}")
            plt.legend(["$v_{c}$", "$V_{c,avg}$"], loc='upper right')
            plt.grid('on')

            # Losses (row 2)
            ax_p = plt.subplot(4, n_caps, 2 * n_caps + i + 1, sharex=ax_i)
            if 'cap' in timeLoss and cid in timeLoss['cap']:
                plt.plot(t_plot[::down2], timeLoss['cap'][cid]['p_L'][::down2])
                plt.plot(t_plot[::down2], np.mean(timeLoss['cap'][cid]['p_L'][::down2]) * np.ones(np.size(t_plot[::down2])))
            plt.ylabel("$p(t)$ (W)")
            plt.title(f"Losses {cid}")
            plt.legend(["$p_{c}$", "$P_{c,avg}$"], loc='upper right')
            plt.grid('on')

            # Temperature (row 3)
            ax_t = plt.subplot(4, n_caps, 3 * n_caps + i + 1, sharex=ax_i)
            if 'cap' in timeTher and cid in timeTher['cap']:
                plt.plot(t_plot[::down2], timeTher['cap'][cid][::down2])
                plt.plot(t_plot[::down2], Ta * np.ones(np.size(t_plot[::down2])))
            plt.ylabel("$\\Theta(t)$ (°C)")
            plt.title(f"Thermal {cid}")
            plt.xlabel(t_label)
            plt.legend(["Cap", "$Ambient$"], loc='upper right')
            plt.grid('on')

    # ==============================================================================
    # 6) Transformer Overview
    # ==============================================================================
    if 'tra' in timeElec and len(timeElec['tra']) > 0:
        plt.figure(figsize=(12, 10))
        plt.suptitle("Transformer Overview: Electrical and Thermal", size=18)
        plt.subplots_adjust(hspace=0.5, wspace=0.3)

        # Electrical: Flux and B
        plt.subplot(4, 1, 1)
        plt.plot(t_plot[::down2], timeElec['tra']['i_m'][::down2])
        plt.title('Magnetizing Current')
        plt.xlabel(t_label)
        plt.ylabel('$i_{m}$ (A)')
        plt.grid('on')

        plt.subplot(4, 1, 2)
        plt.plot(t_plot[::down2], timeElec['tra']['B'][::down2])
        plt.title('Magnetic Flux Density')
        plt.xlabel(t_label)
        plt.ylabel('B (T)')
        plt.grid('on')

        # Losses: Overview
        plt.subplot(4, 1, 3)
        if 'tra' in timeLoss:
            plt.plot(t_plot[::down2], timeLoss['tra']['p_PC'][::down2])
            plt.plot(t_plot[::down2], timeLoss['tra']['p_SC'][::down2])
            plt.plot(t_plot[::down2], timeLoss['tra']['p_CC'][::down2])
            plt.title('Transformer Losses')
            plt.ylabel('Losses (W)')
            plt.legend(["Primary", "Secondary", "Core"], loc='upper right')
            plt.grid('on')

        # Thermal Overview
        plt.subplot(4, 1, 4)
        if 'tra' in timeTher:
            for nid in timeTher['tra']:
                plt.plot(t_plot[::down2], timeTher['tra'][nid][::down2])
        plt.plot(t_plot[::down2], Ta * np.ones(np.size(t_plot[::down2])))
        plt.title('Transformer Temperatures')
        plt.xlabel(t_label)
        plt.ylabel('Temperature (°C)')
        plt.legend(["Primary", "Secondary", "Core", "Ambient"], loc='upper right')
        plt.grid('on')

    plt.show()

    ###################################################################################################################
    # MSG OUT
    ###################################################################################################################
    print("END: Plotting Stationary DAB Waveforms")
