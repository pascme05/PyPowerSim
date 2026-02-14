#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         plotSweep_DAB
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
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import matplotlib.gridspec as gridspec
import matplotlib
matplotlib.use('TkAgg')


#######################################################################################################################
# Function
#######################################################################################################################
def plotSweep_DAB(time, freq, sweep, setup):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("START: Plotting Sweep DAB Waveforms")

    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Init
    # ==============================================================================
    distAc = sweep['Ac']
    distDc = sweep['Dc']
    timeSw = time['Sw']
    timeAc = time['Ac']
    timeDc = time['Dc']
    timeElec = time.get('Elec', {})
    timeLoss = time.get('Loss', {})
    timeTher = time.get('Ther', {})
    freqSw = freq['Sw']
    freqAc = freq['Ac']
    freqDc = freq['Dc']

    # ==============================================================================
    # Parameters
    # ==============================================================================
    fel = setup['Top']['fel']
    fs = setup['Par']['PWM']['fs']
    fsim = setup['Exp']['fsim']
    Q = int(fs / fel)
    Mi = setup['Dat']['stat']['Mi']
    Vdc = setup['Dat']['stat']['Vdc']
    phiDAB = setup['Dat']['stat']['PhiDAB']
    down = int(setup['Dat']['stat']['cyc']) - 2
    down2 = int(fsim / fs / 100)
    if down2 < 1:
        down2 = 1

    # ==============================================================================
    # Variables
    # ==============================================================================
    t = timeSw['t'].values
    f = fsim * np.linspace(0, 0.5, int(len(t) / 2))
    phi_i = np.linspace(0, np.pi / 2, int(setup['Dat']['stat']['W']))

    ###################################################################################################################
    # Pre-Processing
    ###################################################################################################################
    # ==============================================================================
    # Distortion
    # ==============================================================================
    # AC
    # I_ac_pri_thd = thd(timeAc['i_ac_pri'], t, (2 / np.sqrt(2)), 1)
    # I_ac_sec_thd = thd(timeAc['i_ac_sec'], t, (2 / np.sqrt(2)), 1)
    # V_ac_pri_thd = thd(timeAc['v_ac_pri'], t, (2 / np.sqrt(2)), 1)
    # V_ac_sec_thd = thd(timeAc['v_ac_sec'], t, (2 / np.sqrt(2)), 1)

    # DC
    I_dc_pri_thd = thd(timeDc['i_dc_pri'], t, 1, 0)
    I_dc_sec_thd = thd(timeDc['i_dc_sec'], t, 1, 0)
    V_dc_pri_thd = thd(timeDc['v_dc_pri'], t, 1, 0)
    V_dc_sec_thd = thd(timeDc['v_dc_sec'], t, 1, 0)

    # ==============================================================================
    # Time
    # ==============================================================================
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
    # Current
    # ==============================================================================
    plt.figure()
    txt = "Currents DAB for PWM Control with: " \
          + "$V_{dc}$=" + str(Vdc) + "V, " \
          + "$M_{i}$=" + str(Mi) + "$ ,Q$=" + str(Q) \
          + ", Phi=" + f"{math.degrees(phiDAB):.2f}"
    plt.suptitle(txt, size=18)
    plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)

    # ------------------------------------------
    # Phase
    # ------------------------------------------
    # Time
    plt.subplot(2, 3, 1)
    plt.plot(t_plot[::down2], timeAc['i_ac_pri'][::down2])
    plt.plot(t_plot[::down2], timeAc['i_ac_sec'][::down2])
    plt.ylabel("$i_{ac}(t)$ (A)")
    plt.title('Time-domain Currents AC')
    plt.xlabel(t_label)
    pl.legend(["$i_{ac,pri}$", "$i_{ac,sec}$"], loc='upper right')
    plt.grid('on')

    # Frequency
    ax = plt.subplot(2, 3, 2)
    m1, s1, _ = plt.stem(f_plot[::down][0:f_max_plt], freqAc['i_ac_pri'][::down][0:f_max_plt])
    m2, s2, _ = plt.stem(f_plot[::down][0:f_max_plt], freqAc['i_ac_sec'][::down][0:f_max_plt])
    plt.setp(m1, color='C0')
    plt.setp(s1, color='C0')
    plt.setp(m2, color='C1')
    plt.setp(s2, color='C1')
    plt.ylabel("$I_{ac}(f)$ (A)")
    plt.title('Frequency-domain Current AC-Side')
    plt.xlabel(f_label)
    plt.yscale('log')
    plt.ylim(OoM(max(freqAc['i_ac_pri'])) / f_scale, )
    plt.grid('on')

    # Modulation
    plt.subplot(233)
    plt.plot(phi_i, distAc['Pri']['num']['I_a_thd'], label='Primary')
    plt.plot(phi_i, distAc['Sec']['num']['I_a_thd'], label='Secondary')
    if setup['Exp']['plot'] == 2:
        plt.plot(phi_i, distAc['Pri']['ana']['I_a_thd'], 'tab:blue', linestyle="", marker="o", label='Analytical Pri')
        plt.plot(phi_i, distAc['Sec']['ana']['I_a_thd'], 'tab:orange', linestyle="", marker="o", label='Analytical Sec')
    plt.legend(loc='upper right')
    plt.ylim(0, )
    plt.title('Distortion Current AC-Side')
    plt.ylabel("$I_{ac,rms}^{THD}$ (A)")
    plt.xlabel("$Phi$ in (rad)")
    plt.grid('on')

    # ------------------------------------------
    # DC-Link
    # ------------------------------------------
    # Time
    plt.subplot(2, 3, 4)
    plt.plot(t_plot[::down2], timeDc['i_dc_pri'][::down2])
    plt.plot(t_plot[::down2], timeDc['i_dc_sec'][::down2])
    plt.ylabel("$i_{dc}(t)$ (A)")
    plt.title('Time-domain Currents DC-Side')
    plt.xlabel(t_label)
    plt.legend(["$i_{dc,pri}$", "$i_{dc,sec}$", "$I_{dc,pri,avg}$", "$I_{dc,sec,avg}$"], loc='upper right')
    plt.grid('on')

    # Frequency
    ax = plt.subplot(2, 3, 5)
    m1, s1, _ = plt.stem(f_plot[::down][0:f_max_plt], freqDc['i_dc_pri'][::down][0:f_max_plt])
    m2, s2, _ = plt.stem(f_plot[::down][0:f_max_plt], freqDc['i_dc_sec'][::down][0:f_max_plt])
    plt.setp(m1, color='C0')
    plt.setp(s1, color='C0')
    plt.setp(m2, color='C1')
    plt.setp(s2, color='C1')
    plt.ylabel("$I_{dc}(f)$ (A)")
    plt.title('Frequency-domain Current DC-Side')
    plt.xlabel(f_label)
    plt.yscale('log')
    plt.ylim(0.1 / OoM(max(freqDc['i_dc_pri'])), )
    txt = "THD_pri=" + str(round(I_dc_pri_thd * 100, 2)) + "%"
    plt.text(0.70, 0.90, txt, bbox=dict(facecolor='tab:blue', alpha=0.5), transform=ax.transAxes)
    txt = "THD_sec=" + str(round(I_dc_sec_thd * 100, 2)) + "%"
    plt.text(0.70, 0.80, txt, bbox=dict(facecolor='tab:blue', alpha=0.5), transform=ax.transAxes)
    plt.grid('on')

    # Modulation
    plt.subplot(236)
    plt.plot(phi_i, distDc['Pri']['num']['I_dc_thd'], label='Primary')
    plt.plot(phi_i, distDc['Sec']['num']['I_dc_thd'], label='Secondary')
    if setup['Exp']['plot'] == 2:
        plt.plot(phi_i, distDc['Pri']['ana']['I_dc_thd'], 'tab:blue', linestyle="", marker="o", label='Analytical Pri')
        plt.plot(phi_i, distDc['Sec']['ana']['I_dc_thd'], 'tab:orange', linestyle="", marker="o", label='Analytical Sec')
    plt.legend(loc='upper right')
    plt.ylim(0, )
    plt.title('Distortion Current DC-Side')
    plt.ylabel("$I_{dc,rms}^{THD}$ (A)")
    plt.xlabel("$Phi$ in (rad)")
    plt.grid('on')

    # ==============================================================================
    # Voltage
    # ==============================================================================
    plt.figure()
    txt = "Voltages DAB for PWM Control with: " \
          + "$V_{dc}$=" + str(Vdc) + "V, " \
          + "$M_{i}$=" + str(Mi) + "$ ,Q$=" + str(Q) \
          + ", Phi=" + f"{math.degrees(phiDAB):.2f}"
    plt.suptitle(txt, size=18)
    plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)

    # ------------------------------------------
    # Phase
    # ------------------------------------------
    # Time
    plt.subplot(2, 3, 1)
    plt.plot(t_plot[::down2], timeAc['v_ac_pri'][::down2])
    plt.plot(t_plot[::down2], timeAc['v_L'][::down2])
    plt.plot(t_plot[::down2], timeAc['v_ac_sec_ref'][::down2])
    plt.ylabel("$v_{ac}(t)$ (V)")
    plt.title('Time-domain Voltages AC-Side')
    plt.xlabel(t_label)
    plt.legend(["$v_{ac,pri}$", "$v_{ac,Lk}$", "$v_{ac,sec,ref}$"], loc='upper right')
    plt.grid('on')

    # Frequency
    ax = plt.subplot(2, 3, 2)
    m1, s1, _ = plt.stem(f_plot[::down][0:f_max_plt], freqAc['v_ac_pri'][::down][0:f_max_plt])
    m2, s2, _ = plt.stem(f_plot[::down][0:f_max_plt], freqAc['v_L'][::down][0:f_max_plt])
    plt.setp(m1, color='C0')
    plt.setp(s1, color='C0')
    plt.setp(m2, color='C1')
    plt.setp(s2, color='C1')
    plt.ylabel("$V_{ac}(f)$ (V)")
    plt.title('Frequency-domain Voltages AC-Side')
    plt.xlabel(f_label)
    plt.yscale('log')
    plt.ylim(OoM(max(freqAc['v_ac_pri'])) / f_scale, )
    plt.grid('on')

    # Modulation
    plt.subplot(233)
    plt.plot(phi_i, distAc['Pri']['num']['V_a_thd'], label='Primary')
    plt.plot(phi_i, distAc['Sec']['num']['V_a_thd'], label='Secondary')
    if setup['Exp']['plot'] == 2:
        plt.plot(phi_i, distAc['Pri']['ana']['V_a_thd'], 'tab:blue', linestyle="", marker="o", label='Analytical Pri')
        plt.plot(phi_i, distAc['Sec']['ana']['V_a_thd'], 'tab:orange', linestyle="", marker="o", label='Analytical Sec')
    plt.legend(loc='upper right')
    plt.title('Distortion Voltage AC-Side')
    plt.ylabel("$V_{a,rms}^{THD}$ (V)")
    plt.xlabel("$Phi$ in (rad)")
    plt.grid('on')

    # ------------------------------------------
    # DC-Link
    # ------------------------------------------
    # Time
    plt.subplot(2, 3, 4)
    plt.plot(t_plot[::down2], timeDc['v_dc_pri'][::down2])
    plt.plot(t_plot[::down2], timeDc['v_dc_sec'][::down2])
    plt.ylabel("$v_{dc}(t)$ (V)")
    plt.title('Time-domain Voltages DC-Side')
    plt.xlabel(t_label)
    plt.legend(["$v_{dc,pri}$", "$v_{dc,sec}$"], loc='upper right')
    plt.grid('on')

    # Frequency
    ax = plt.subplot(2, 3, 5)
    m1, s1, _ = plt.stem(f_plot[::down][0:f_max_plt], freqDc['v_dc_pri'][::down][0:f_max_plt])
    m2, s2, _ = plt.stem(f_plot[::down][0:f_max_plt], freqDc['v_dc_sec'][::down][0:f_max_plt])
    plt.setp(m1, color='C0')
    plt.setp(s1, color='C0')
    plt.setp(m2, color='C1')
    plt.setp(s2, color='C1')
    plt.ylabel("$V_{dc}(f)$ (V)")
    plt.title('Frequency-domain Voltages DC-Side')
    plt.xlabel(f_label)
    plt.yscale('log')
    plt.ylim(0.1 / OoM(max(freqDc['v_dc_pri'])), )
    txt = "THD_pri=" + str(round(V_dc_pri_thd * 100, 2)) + "%"
    plt.text(0.7, 0.90, txt, bbox=dict(facecolor='tab:blue', alpha=0.5), transform=ax.transAxes)
    txt = "THD_sec=" + str(round(V_dc_sec_thd * 100, 2)) + "%"
    plt.text(0.7, 0.80, txt, bbox=dict(facecolor='tab:blue', alpha=0.5), transform=ax.transAxes)
    plt.grid('on')

    # Modulation
    plt.subplot(236)
    plt.plot(phi_i, distDc['Pri']['num']['V_dc_thd'], label='Primary')
    plt.plot(phi_i, distDc['Sec']['num']['V_dc_thd'], label='Secondary')
    if setup['Exp']['plot'] == 2:
        plt.plot(phi_i, distDc['Pri']['ana']['V_dc_thd'], 'tab:blue', linestyle="", marker="o", label='Analytical Pri')
        plt.plot(phi_i, distDc['Sec']['ana']['V_dc_thd'], 'tab:orange', linestyle="", marker="o", label='Analytical Sec')
    plt.legend(loc='upper right')
    plt.title('Distortion Voltage DC-Side')
    plt.ylabel("$V_{dc,rms}^{THD}$ (V)")
    plt.xlabel("$Phi$ in (rad)")
    plt.grid('on')
    plt.show()

    ###################################################################################################################
    # MSG OUT
    ###################################################################################################################
    print("END: Plotting DAB Sweep Results")
