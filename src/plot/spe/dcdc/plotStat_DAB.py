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
    down2 = int(fsim / fs / 200)
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
    f_max = f[-1]
    if f_max < 1e3:
        f_plot = f
        f_label = "$f$ in (Hz)"
    elif f_max < 1e6:
        f_plot = f / 1e3
        f_label = "$f$ in (kHz)"
    else:
        f_plot = f / 1e6
        f_label = "$f$ in (MHz)"

    ###################################################################################################################
    # 1) Switching strategies
    ###################################################################################################################
    gs = gridspec.GridSpec(3, 2)
    pl.figure(figsize=(12, 10))
    txt = "DAB Switching Strategies: " + "$M_{i}$=" + str(Mi) + "$ ,Q$=" + str(Q) + ", Ntr=" + str(n_tr) + ", Phi=" + str(phiDAB)
    pl.suptitle(txt, size=18)
    pl.subplots_adjust(hspace=0.4, wspace=0.3, left=0.1, right=0.9, top=0.9, bottom=0.1)

    # ------------------------------------------
    # Modulation
    # ------------------------------------------
    ax1 = pl.subplot(gs[0, :])
    ax2 = ax1.twinx()
    ax1.plot(t_plot[::down2], timeSw['cA'][::down2], 'tab:blue', label='$c_{pri}$')
    ax1.plot(t_plot[::down2], timeSw['cB'][::down2], 'tab:orange', label='$c_{sec}$')
    ax2.plot(t_plot[::down2], timeSw['v_a_ref'][::down2] / (Vdc / 2), color='tab:blue', linestyle='--', label='$v_{pri}^{*}$')
    ax2.plot(t_plot[::down2], timeSw['v_b_ref'][::down2] / (Vdc / 2), color='tab:orange', linestyle='--', label='$v_{sec}^{*}$')
    pl.title('Carrier and Reference Waveforms')
    pl.xlabel(t_label)
    ax1.set_ylabel('c(t)')
    ax2.set_ylabel('v(t) (p.u)')
    ax1.legend(loc='upper left')
    ax2.legend(loc='upper right')
    pl.grid('on')

    # ------------------------------------------
    # Switching Waveforms
    # ------------------------------------------
    # Primary
    pl.subplot(gs[1, 0])
    pl.plot(t_plot[::down2], timeSw['sA'][::down2], label='$s_{A}$')
    pl.plot(t_plot[::down2], timeSw['sB'][::down2], label='$s_{B}$')
    pl.title('Primary Bridge Switching')
    pl.xlabel(t_label)
    pl.ylabel("$s(t)$ (p.u)")
    pl.legend(loc='upper right')
    pl.grid('on')

    # Secondary
    pl.subplot(gs[1, 1])
    pl.plot(t_plot[::down2], timeSw['sC'][::down2], label='$s_{C}$')
    pl.plot(t_plot[::down2], timeSw['sD'][::down2], label='$s_{D}$')
    pl.title('Secondary Bridge Switching')
    pl.xlabel(t_label)
    pl.ylabel("$s(t)$ (p.u)")
    pl.legend(loc='upper right')
    pl.grid('on')

    # ------------------------------------------
    # Spectrum
    # ------------------------------------------
    # Primary
    pl.subplot(gs[2, 0])
    pl.stem(f_plot[::down][0:50], freqSw['Sa'][::down][0:50])
    pl.xlim(0, 50)
    pl.title('Primary Switching Spectrum')
    pl.xlabel(f_label)
    pl.ylabel("$S_{p}(f)$ (p.u)")
    pl.yscale('log')
    pl.ylim(0.0001, )
    pl.grid('on')

    # Secondary
    pl.subplot(gs[2, 1])
    if 'Sb' in freqSw:
        pl.stem(f_plot[::down][0:50], freqSw['Sb'][::down][0:50])
    elif 'Sc' in freqSw:
        pl.stem(f_plot[::down][0:50], freqSw['Sc'][::down][0:50])
    pl.xlim(0, 50)
    pl.title('Secondary Switching Spectrum')
    pl.xlabel(f_label)
    pl.ylabel("$S_{s}(f)$ (p.u)")
    pl.yscale('log')
    pl.ylim(0.0001, )
    pl.grid('on')

    ###################################################################################################################
    # 2) Primary Side Current and Voltage
    ###################################################################################################################
    plt.figure(figsize=(12, 8))
    plt.suptitle("DAB Primary Side: Current and Voltage", size=18)
    plt.subplots_adjust(hspace=0.4, wspace=0.3)

    # Voltage Time AC
    plt.subplot(2, 2, 1)
    plt.plot(t_plot[::down2], timeAc['v_ac_pri'][::down2])
    plt.title('Primary Voltage (Time Domain)')
    plt.xlabel(t_label)
    plt.ylabel("$v_{p}(t)$ (V)")
    plt.grid('on')

    # Voltage Freq
    ax = plt.subplot(2, 2, 2)
    v_p_spec = freqAc.get('v_ac_pri', freqAc.get('V_a'))
    plt.stem(f_plot[::down][0:50], v_p_spec[::down][0:50])
    plt.title('Primary Voltage (Frequency Domain)')
    plt.xlabel("$f/f_{1}$")
    plt.ylabel("$V_{p}(f)$ (V)")
    plt.yscale('log')
    pl.ylim(0.0001, )
    plt.grid('on')

    # Current Time AC
    plt.subplot(2, 2, 3)
    plt.plot(t_plot[::down2], timeAc['i_ac_pri'][::down2])
    plt.title('Primary Current (Time Domain)')
    plt.xlabel(t_label)
    plt.ylabel("$i_{p}(t)$ (A)")
    plt.grid('on')

    # Current Freq
    ax = plt.subplot(2, 2, 4)
    i_p_spec = freqAc.get('i_ac_pri', freqAc.get('I_a'))
    plt.stem(f_plot[::down][0:50], i_p_spec[::down][0:50])
    plt.title('Primary Current (Frequency Domain)')
    plt.xlabel("$f/f_{1}$")
    plt.ylabel("$I_{p}(f)$ (A)")
    plt.yscale('log')
    plt.grid('on')

    ###################################################################################################################
    # 3) Secondary Side Current and Voltage
    ###################################################################################################################
    plt.figure(figsize=(12, 8))
    plt.suptitle("DAB Secondary Side: Current and Voltage", size=18)
    plt.subplots_adjust(hspace=0.4, wspace=0.3)

    # Voltage Time
    plt.subplot(2, 2, 1)
    plt.plot(t_plot[::down2], timeAc['v_ac_sec'][::down2])
    plt.title('Secondary Voltage (Time Domain)')
    plt.xlabel(t_label)
    plt.ylabel("$v_{s}(t)$ (V)")
    plt.grid('on')

    # Voltage Freq
    ax = plt.subplot(2, 2, 2)
    v_s_spec = freqAc.get('v_ac_sec', freqAc.get('V_b'))
    if v_s_spec is None and 'v_b0' in freqAc: v_s_spec = freqAc['v_b0']
    if v_s_spec is not None:
        plt.stem(f_plot[::down][0:50], v_s_spec[::down][0:50])
    plt.title('Secondary Voltage (Frequency Domain)')
    plt.xlabel("$f/f_{1}$")
    plt.ylabel("$V_{s}(f)$ (V)")
    plt.yscale('log')
    plt.grid('on')

    # Current Time
    plt.subplot(2, 2, 3)
    plt.plot(t_plot[::down2], timeAc['i_ac_sec'][::down2])
    plt.title('Secondary Current (Time Domain)')
    plt.xlabel(t_label)
    plt.ylabel("$i_{s}(t)$ (A)")
    plt.grid('on')

    # Current Freq
    ax = plt.subplot(2, 2, 4)
    i_s_spec = freqAc.get('i_ac_sec', freqAc.get('I_b'))
    if i_s_spec is not None:
        plt.stem(f_plot[::down][0:50], i_s_spec[::down][0:50])
    plt.title('Secondary Current (Frequency Domain)')
    plt.xlabel("$f/f_{1}$")
    plt.ylabel("$I_{s}(f)$ (A)")
    plt.yscale('log')
    plt.grid('on')

    ###################################################################################################################
    # 4) Overview power losses switches and capacitor
    ###################################################################################################################
    if 'sw' in timeLoss:
        plt.figure(figsize=(12, 8))
        plt.suptitle("Power Losses Overview: Switches and Capacitors", size=18)
        plt.subplots_adjust(hspace=0.4, wspace=0.3)

        # Switches Losses
        plt.subplot(2, 2, 1)
        p_sw_pri = []
        p_sw_sec = []
        labels_pri = []
        labels_sec = []
        for i in range(1, 5):
            sid = 'S' + str(i)
            if sid in timeLoss['sw']:
                p_sw_pri.append(np.mean(timeLoss['sw'][sid]['p_L']))
                labels_pri.append(sid)
        for i in range(5, 9):
            sid = 'S' + str(i)
            if sid in timeLoss['sw']:
                p_sw_sec.append(np.mean(timeLoss['sw'][sid]['p_L']))
                labels_sec.append(sid)

        plt.bar(labels_pri + labels_sec, p_sw_pri + p_sw_sec)
        plt.title('Average Switch Losses')
        plt.ylabel('Losses (W)')
        plt.grid('on', axis='y')

        # Capacitor Losses
        plt.subplot(2, 2, 2)
        p_cap = []
        labels_cap = []
        if 'cap' in timeLoss:
            for cid in timeLoss['cap']:
                p_cap.append(np.mean(timeLoss['cap'][cid]['p_L']))
                labels_cap.append(cid)
        plt.bar(labels_cap, p_cap, color='tab:orange')
        plt.title('Average Capacitor Losses')
        plt.ylabel('Losses (W)')
        plt.grid('on', axis='y')

        # Losses over time (Summed)
        plt.subplot(2, 1, 2)
        p_total_sw = np.zeros(len(t))
        for sid in timeLoss['sw']:
            p_total_sw += timeLoss['sw'][sid]['p_L'].values
        plt.plot(t_plot[::down2], p_total_sw[::down2], label='Total Switch Losses')
        if 'cap' in timeLoss:
            p_total_cap = np.zeros(len(t))
            for cid in timeLoss['cap']:
                p_total_cap += timeLoss['cap'][cid]['p_L'].values
            plt.plot(t_plot[::down2], p_total_cap[::down2], label='Total Cap Losses')
        plt.title('Instantaneous Losses')
        plt.xlabel(t_label)
        plt.ylabel('Power (W)')
        plt.legend()
        plt.grid('on')

    ###################################################################################################################
    # 5) Thermal switches and capacitor
    ###################################################################################################################
    if 'sw' in timeTher:
        plt.figure(figsize=(12, 8))
        plt.suptitle("Thermal Overview: Switches and Capacitors", size=18)
        plt.subplots_adjust(hspace=0.4, wspace=0.3)

        # Primary Switches
        plt.subplot(2, 2, 1)
        for i in range(1, 5):
            sid = 'T' + str(i)
            if sid in timeTher['sw']:
                plt.plot(t_plot[::down2], timeTher['sw'][sid][::down2], label=sid)
        plt.title('Primary Switches Temperature')
        plt.xlabel(t_label)
        plt.ylabel('Temperature (°C)')
        plt.legend()
        plt.grid('on')

        # Secondary Switches
        plt.subplot(2, 2, 2)
        for i in range(5, 9):
            sid = 'T' + str(i)
            if sid in timeTher['sw']:
                plt.plot(t_plot[::down2], timeTher['sw'][sid][::down2], label=sid)
        plt.title('Secondary Switches Temperature')
        plt.xlabel(t_label)
        plt.ylabel('Temperature (°C)')
        plt.legend()
        plt.grid('on')

        # Capacitor Thermal
        plt.subplot(2, 2, 3)
        if 'cap' in timeTher:
            for cid in timeTher['cap']:
                plt.plot(t_plot[::down2], timeTher['cap'][cid][::down2], label=cid)
        plt.title('Capacitor Temperature')
        plt.xlabel(t_label)
        plt.ylabel('Temperature (°C)')
        plt.legend()
        plt.grid('on')

        # Heatsink/Case (if available)
        plt.subplot(2, 2, 4)
        for i in range(1, 9, 4): # S1 and S5
            sid = 'C' + str(i)
            if sid in timeTher['sw']:
                plt.plot(t_plot[::down2], timeTher['sw'][sid][::down2], label='Case '+str(i))
        plt.title('Case/Heatsink Temperature')
        plt.xlabel(t_label)
        plt.ylabel('Temperature (°C)')
        plt.legend()
        plt.grid('on')

    ###################################################################################################################
    # 6) Overview Trafo electrical and thermal
    ###################################################################################################################
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
