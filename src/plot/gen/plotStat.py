#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         plotStat
# Date:         11.05.2024
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
def plotStat(time, freq, setup):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("START: Plotting Stationary Waveforms")

    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # id
    # ==============================================================================
    if 'Elec' in time and 'sw' in time['Elec']:
        def _key(x):
            try:
                return int(x[1:])
            except:
                return 0
        id1 = sorted(list(time['Elec']['sw'].keys()), key=_key)
    else:
        id1 = ['S1', 'S2', 'S3', 'S4', 'S5', 'S6']

    # ==============================================================================
    # Init
    # ==============================================================================
    # ------------------------------------------
    # Time
    # ------------------------------------------
    timeSw = time['Sw']
    timeAc = time['Ac']
    timeDc = time['Dc']
    timeElec = time['Elec']
    timeLoss = time['Loss']
    timeTher = time['Ther']

    # ------------------------------------------
    # Frequency
    # ------------------------------------------
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
    R = setup['Top']['R']
    L = setup['Top']['L']
    Mi = setup['Dat']['stat']['Mi']
    Vdc = setup['Dat']['stat']['Vdc']
    phiE = setup['Top']['phiE']
    down = int(setup['Dat']['stat']['cyc']) - 2
    down2 = int(fsim / fs / 200)
    Ta = setup['Dat']['stat']['Tc']
    if down2 < 1:
        down2 = 1

    # ==============================================================================
    # Variables
    # ==============================================================================
    t = timeSw['t'].values
    f = fsim * np.linspace(0, 0.5, int(len(t) / 2)) / fel

    ###################################################################################################################
    # Pre-Processing
    ###################################################################################################################
    # ==============================================================================
    # General
    # ==============================================================================
    # ------------------------------------------
    # Limits
    # ------------------------------------------
    K = int(np.round((t[-1] - t[0]) * fel))
    start = int((len(t) - 1) / K) * (K - 1)
    ende = len(t)

    # ------------------------------------------
    # Change time
    # ------------------------------------------
    t = t[start:ende]
    timeSw = timeSw[:][start:ende]
    for c1 in timeAc:
        timeAc[c1] = timeAc[c1][start:ende]
    for c1 in timeDc:
        timeDc[c1] = timeDc[c1][start:ende]
    for c1 in timeElec:
        for c2 in timeElec[c1]:
            timeElec[c1][c2] = timeElec[c1][c2][start:ende]
    for c1 in timeLoss:
        for c2 in timeLoss[c1]:
            timeLoss[c1][c2] = timeLoss[c1][c2][start:ende]
    for c1 in timeTher:
        for c2 in timeTher[c1]:
            timeTher[c1][c2] = timeTher[c1][c2][start:ende]

    # ==============================================================================
    # Load angle Total
    # ==============================================================================
    Y = fft(timeAc['v_a'])
    phiV = np.angle(Y)[1]
    Y = fft(timeAc['i_a'])
    phiI = np.angle(Y)[1]
    phi = phiV - phiI + 2 * np.pi
    while phi > 2 * np.pi:
        phi = phi - 2 * np.pi

    # ==============================================================================
    # Load angle RL
    # ==============================================================================
    angZ = math.atan2(2 * np.pi * fel * L, R)

    # ==============================================================================
    # Distortion
    # ==============================================================================
    I_ph_thd = thd(timeAc['i_a'], t, (2 / np.sqrt(2)), 1)
    V_ph_thd = thd(timeAc['v_a'], t, (2 / np.sqrt(2)), 1)
    I_dc_thd = thd(timeDc['i_dc'], t, 1, 0)
    V_dc_thd = thd(timeDc['v_dc'], t, 1, 0)

    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    # ==============================================================================
    # Switching Function
    # ==============================================================================
    gs = gridspec.GridSpec(3, 2)
    pl.figure()
    txt = "Modulation Functions for: " + "$M_{i}$=" + str(Mi) + "$ ,Q$=" + str(Q) + ", Sampling: " + str(
        setup['Par']['PWM']['samp']) + ", Update: " + str(setup['Par']['PWM']['upd']) + " and Edge Trigger: " + str(
        setup['Par']['PWM']['tri'])
    pl.suptitle(txt, size=18)
    pl.subplots_adjust(hspace=0.35, wspace=0.20, left=0.075, right=0.925, top=0.90, bottom=0.075)

    # ------------------------------------------
    # Modulation
    # ------------------------------------------
    ax1 = pl.subplot(gs[0, :])
    ax2 = ax1.twinx()
    ax1.plot(t[::down2], timeSw['c'][::down2], 'black')
    ax2.plot(t[::down2], timeSw['v_a_ref'][::down2] / (Vdc / 2), color='tab:blue')
    ax2.plot(t[::down2], timeSw['xA'][::down2], '--', color='tab:blue')
    ax2.plot(t[::down2], timeSw['n0'][::down2], '-.', color='tab:blue')
    pl.title('Carrier and Reference Waveform')
    pl.xlabel('t in (sec)')
    ax1.set_ylabel('c(t)/s(t)', color='black')
    ax2.set_ylabel('v(t) (p.u)', color='tab:blue')
    pl.legend(["$c_{a}$", "$v_{a}^{*}$", "$v_{a0}^{*}$", "$v_{n0}^{*}$"], loc='upper right')
    pl.grid('on')

    # ------------------------------------------
    # Switching Waveform
    # ------------------------------------------
    # Time-Domain
    pl.subplot(gs[1, 0])
    pl.plot(t[::down2], timeSw['sA'][::down2])
    pl.title('Time-domain Switching Functions')
    pl.xlabel('t in (sec)')
    pl.ylabel("$s_{a}(t)$ (p.u)")
    pl.legend(["$s_{a}$"], loc='upper right')
    pl.grid('on')

    # Freq-Domain
    pl.subplot(gs[2, 0])
    pl.stem(f[::down][0:50], freqSw['Sa'][::down][0:50])
    pl.xlim(0, 50)
    pl.title('Frequency-domain Switching Function')
    pl.xlabel("$f/f_{1}$ (Hz/Hz)")
    pl.ylabel("$S_{a}(f)$ (p.u)")
    pl.yscale('log')
    pl.ylim(0.0001, )
    pl.grid('on')

    # ------------------------------------------
    # Reference Waveform
    # ------------------------------------------
    # Time-Domain
    pl.subplot(gs[1, 1])
    pl.plot(t[::down2], timeSw['xAs'][::down2])
    pl.title('Time-domain Sampled References')
    pl.xlabel('t in (sec)')
    pl.ylabel("$x_{a}(t)$ (p.u)")
    pl.legend(["$x_{a}$"], loc='upper right')
    pl.grid('on')

    # Freq-Domain
    pl.subplot(gs[2, 1])
    pl.stem(f[::down][0:50], freqSw['Xas'][::down][0:50])
    pl.xlim(0, 50)
    pl.title('Frequency-domain Sampled Reference')
    pl.xlabel("$f/f_{1}$ (Hz/Hz)")
    pl.ylabel("$X_{a}(f)$ (p.u)")
    pl.yscale('log')
    pl.ylim(0.0001, )
    pl.grid('on')

    # ==============================================================================
    # Current
    # ==============================================================================
    plt.figure()
    txt = "Currents for PWM Control with: " + "$V_{dc}$=" + str(Vdc) + "V, " + "$M_{i}$=" + str(
        Mi) + "$ ,Q$=" + str(Q) + ", $\phi_{RL}=$" + str(int(math.degrees(angZ))) + "deg" + ", $\phi_{E}=$" + str(
        int(phiE)) + "deg" + ", $\phi_{VI}=$" + str(int(math.degrees(phi))) + "deg"
    plt.suptitle(txt, size=18)
    plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)

    # ------------------------------------------
    # Phase
    # ------------------------------------------
    # Time
    plt.subplot(2, 2, 1)
    plt.plot(t[::down2], timeAc['i_a'][::down2])
    plt.ylabel("$i_{a}(t)$ (A)")
    plt.title('Time-domain Currents AC-Side')
    plt.xlabel('time in (sec)')
    pl.legend(["$i_{a}$"], loc='upper right')
    plt.grid('on')

    # Frequency
    ax = plt.subplot(2, 2, 2)
    plt.stem(f[::down][0:50], freqAc['I_a'][::down][0:50])
    plt.ylabel("$I_{a}(f)$ (A)")
    plt.xlim(0, 50)
    plt.title('Frequency-domain Current AC-Side')
    plt.xlabel("$f/f_{1}$ (Hz/Hz)")
    plt.yscale('log')
    plt.ylim(0.1 / OoM(max(freqAc['I_a'])), )
    txt = "THD=" + str(round(I_ph_thd * 100, 2)) + "%"
    plt.text(0.75, 0.90, txt, bbox=dict(facecolor='tab:blue', alpha=0.5), transform=ax.transAxes)
    plt.grid('on')

    # ------------------------------------------
    # DC-Link
    # ------------------------------------------
    # Time
    plt.subplot(2, 2, 3)
    plt.plot(t[::down2], timeDc['i_dc'][::down2], t[::down2],
             np.mean(timeDc['i_dc']) * np.ones(np.size(timeDc['i_dc'][::down2])), '--')
    plt.ylabel("$i_{dc}(t)$ (A)")
    plt.title('Time-domain Currents DC-Side')
    plt.xlabel('time in (sec)')
    plt.legend(["$i_{dc}$", "$I_{dc,avg}$"], loc='upper right')
    plt.grid('on')

    # Frequency
    ax = plt.subplot(2, 2, 4)
    plt.stem(f[::down][0:50], freqDc['I_dc'][::down][0:50])
    plt.ylabel("$I_{dc}(f)$ (A)")
    plt.xlim(0, 50)
    plt.title('Frequency-domain Current DC-Side')
    plt.xlabel("$f/f_{1}$ (Hz/Hz)")
    plt.yscale('log')
    plt.ylim(0.1 / OoM(max(freqDc['I_dc'])), )
    txt = "THD=" + str(round(I_dc_thd * 100, 2)) + "%"
    plt.text(0.75, 0.90, txt, bbox=dict(facecolor='tab:blue', alpha=0.5), transform=ax.transAxes)
    plt.grid('on')

    # ==============================================================================
    # Voltage
    # ==============================================================================
    plt.figure()
    txt = "Voltages for PWM Control with: " + "$V_{dc}$=" + str(Vdc) + "V, " + "$M_{i}$=" + str(
        Mi) + ", $\phi_{RL}=$" + str(int(math.degrees(angZ))) + "deg" + ", $\phi_{E}=$" + str(
        int(phiE)) + "deg" + ", $\phi_{VI}=$" + str(int(math.degrees(phi))) + "deg"
    plt.suptitle(txt, size=18)
    plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)

    # ------------------------------------------
    # Phase
    # ------------------------------------------
    # Time
    plt.subplot(2, 2, 1)
    plt.plot(t[::down2], timeAc['v_a'][::down2], t[::down2], timeAc['v_a0'][::down2], t[::down2],
             timeAc['v_a_out'][::down2], t[::down2], timeSw['e_a'][::down2])
    plt.ylabel("$v_{a}(t)$ (V)")
    plt.title('Time-domain Voltages AC-Side')
    plt.xlabel('time in (sec)')
    plt.legend(["$v_{a}(t)$", "$v_{a0}(t)$", "$v_{a,out}(t)$", "$e_{a}(t)$"], loc='upper right')
    plt.grid('on')

    # Frequency
    ax = plt.subplot(2, 2, 2)
    plt.stem(f[::down][0:50], freqAc['V_a'][::down][0:50])
    plt.ylabel("$V_{a}(f)$ (V)")
    plt.xlim(0, 50)
    plt.title('Frequency-domain Voltages AC-Side')
    plt.xlabel("$f/f_{1}$ (Hz/Hz)")
    plt.yscale('log')
    plt.ylim(0.1 / OoM(max(freqAc['V_a'])), )
    txt = "THD=" + str(round(V_ph_thd * 100, 2)) + "%"
    plt.text(0.75, 0.90, txt, bbox=dict(facecolor='tab:blue', alpha=0.5), transform=ax.transAxes)
    plt.grid('on')

    # ------------------------------------------
    # DC-Link
    # ------------------------------------------
    # Time
    plt.subplot(2, 2, 3)
    plt.plot(t[::down2], timeDc['v_in'][::down2], t[::down2], timeDc['v_dc'][::down2], t[::down2],
             np.mean(timeDc['v_dc']) * np.ones(np.size(timeDc['v_dc'][::down2])), '--')
    plt.ylabel("$v_{dc}(t)$ (V)")
    plt.title('Time-domain Voltages DC-Side')
    plt.xlabel('time in (sec)')
    plt.legend(["$v_{in}$", "$v_{dc}$", "$V_{dc,avg}$"], loc='upper right')
    plt.grid('on')

    # Frequency
    ax = plt.subplot(2, 2, 4)
    plt.stem(f[::down][0:50], freqDc['V_dc'][::down][0:50])
    plt.ylabel("$V_{dc}(f)$ (A)")
    plt.xlim(0, 50)
    plt.title('Frequency-domain Voltages DC-Side')
    plt.xlabel("$f/f_{1}$ (Hz/Hz)")
    plt.yscale('log')
    plt.ylim(0.1 / OoM(max(freqDc['V_dc'])), )
    txt = "THD=" + str(round(V_dc_thd * 100, 2)) + "%"
    plt.text(0.75, 0.90, txt, bbox=dict(facecolor='tab:blue', alpha=0.5), transform=ax.transAxes)
    plt.grid('on')

    # ==============================================================================
    # Time-domain Switches
    # ==============================================================================
    plt.figure()
    txt = "Time domain switches for PWM control with: " + "$V_{dc}$=" + str(Vdc) + "V, " + "$M_{i}$=" + str(
        Mi) + ", $\phi_{RL}=$" + str(int(math.degrees(angZ))) + "deg" + ", $\phi_{E}=$" + str(
        int(phiE)) + "deg" + ", $\phi_{VI}=$" + str(int(math.degrees(phi))) + "deg"
    plt.suptitle(txt, size=18)
    plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)

    # ------------------------------------------
    # Losses
    # ------------------------------------------
    for i in range(0, 2):
        plt.subplot(4, 2, 1 + i)
        plt.plot(t[::down2], timeElec['sw'][id1[i]]['i_T'][::down2], t[::down2], timeElec['sw'][id1[i]]['i_D'][::down2])
        plt.ylabel("$i(t)$ (A)")
        txt = 'Time-domain Currents Switch ' + str(id1[i])
        plt.title(txt)
        plt.xticks([], [])
        plt.legend(["$i_{T}$", "$i_{D}$"], loc='upper right')
        plt.grid('on')

        # Voltage
        plt.subplot(4, 2, 3 + i)
        plt.plot(t[::down2], timeElec['sw'][id1[i]]['v_T'][::down2], t[::down2], timeElec['sw'][id1[i]]['v_D'][::down2])
        plt.ylabel("$v(t)$ (V)")
        txt = 'Time-domain Voltages Switch ' + str(id1[i])
        plt.title(txt)
        plt.xticks([], [])
        plt.legend(["$v_{T}$", "$v_{D}$"], loc='upper right')
        plt.grid('on')

        # Losses
        plt.subplot(4, 2, 5 + i)
        plt.plot(t[::down2], timeLoss['sw'][id1[i]]['p_T_c'][::down2], t[::down2], timeLoss['sw'][id1[i]]['p_D_c'][::down2])
        plt.ylabel("$p(t)$ (W)")
        txt = 'Time-domain Conduction Losses Switch ' + str(id1[i])
        plt.title(txt)
        plt.xticks([], [])
        plt.legend(["$p_{T}$", "$p_{D}$"], loc='upper right')
        plt.grid('on')

        # Losses
        plt.subplot(4, 2, 7 + i)
        plt.plot(t[::down2], timeLoss['sw'][id1[i]]['p_T_s'][::down2], t[::down2], timeLoss['sw'][id1[i]]['p_D_s'][::down2])
        plt.ylabel("$p(t)$ (W)")
        txt = 'Time-domain Switching Losses Switch ' + str(id1[i])
        plt.title(txt)
        plt.xlabel('time in (sec)')
        plt.legend(["$p_{T}$", "$p_{D}$"], loc='upper right')
        plt.grid('on')

    # ------------------------------------------
    # Thermal
    # ------------------------------------------
    plt.figure()
    plt.subplot(2, 2, 1)
    plt.plot(t[::down2], timeLoss['sw']['S1']['p_T'][::down2], t[::down2], timeLoss['sw']['S1']['p_D'][::down2])
    plt.ylabel("$p(t)$ (W)")
    txt = 'Time-domain Total Losses Switch S1'
    plt.title(txt)
    plt.xticks([], [])
    plt.legend(["$p_{T}$", "$p_{D}$"], loc='upper right')
    plt.grid('on')

    # Temperature
    plt.subplot(2, 2, 3)
    plt.plot(t[::down2], timeTher['sw']['T1'][::down2], t[::down2], timeTher['sw']['D1'][::down2], t[::down2], timeTher['sw']['C1'][::down2], t[::down2],
             Ta * np.ones(np.size(t[::down2])))
    plt.ylabel("$\Theta(t)$ (°C)")
    txt = 'Time-domain Thermal Switch S1'
    plt.title(txt)
    plt.xlabel('time in (sec)')
    plt.legend(["$\Theta_{T}$", "$\Theta_{D}$", "$\Theta_{C}$", "$\Theta_{A}$"], loc='upper right')
    plt.grid('on')

    # Losses
    plt.subplot(2, 2, 2)
    plt.plot(t[::down2], timeLoss['sw']['S2']['p_T'][::down2], t[::down2], timeLoss['sw']['S2']['p_D'][::down2])
    plt.ylabel("$p(t)$ (W)")
    txt = 'Time-domain Total Losses Switch S2'
    plt.title(txt)
    plt.xticks([], [])
    plt.legend(["$p_{T}$", "$p_{D}$"], loc='upper right')
    plt.grid('on')

    # Temperature
    plt.subplot(2, 2, 4)
    plt.plot(t[::down2], timeTher['sw']['T2'][::down2], t[::down2], timeTher['sw']['D2'][::down2], t[::down2],
             timeTher['sw']['C2'][::down2], t[::down2],
             Ta * np.ones(np.size(t[::down2])))
    plt.ylabel("$\Theta(t)$ (°C)")
    txt = 'Time-domain Thermal Switch S2'
    plt.title(txt)
    plt.xlabel('time in (sec)')
    plt.legend(["$\Theta_{T}$", "$\Theta_{D}$", "$\Theta_{C}$", "$\Theta_{A}$"], loc='upper right')
    plt.grid('on')

    # ==============================================================================
    # Time-domain Capacitor
    # ==============================================================================
    plt.figure()
    txt = "Time domain capacitor for PWM control with: " + "$V_{dc}$=" + str(Vdc) + "V, " + "$M_{i}$=" + str(
        Mi) + ", $\phi_{RL}=$" + str(int(math.degrees(angZ))) + "deg" + ", $\phi_{E}=$" + str(
        int(phiE)) + "deg" + ", $\phi_{VI}=$" + str(int(math.degrees(phi))) + "deg"
    plt.suptitle(txt, size=18)
    plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)

    # ------------------------------------------
    # DC Link
    # ------------------------------------------
    # Current
    plt.subplot(4, 1, 1)
    plt.plot(t[::down2], timeDc['i_c'][::down2])
    plt.ylabel("$i(t)$ (A)")
    plt.title('Time-domain Currents DC-Link Capacitor')
    plt.grid('on')

    # Voltage
    plt.subplot(4, 1, 2)
    plt.plot(t[::down2], timeDc['v_dc'][::down2])
    plt.ylabel("$v(t)$ (V)")
    plt.title('Time-domain Voltages DC-Link Capacitor')
    plt.grid('on')

    # Losses
    plt.subplot(4, 1, 3)
    plt.plot(t[::down2], timeLoss['cap']['C1']['p_L'][::down2])
    plt.ylabel("$p(t)$ (W)")
    plt.title('Time-domain Losses DC-Link Capacitor')
    plt.grid('on')

    # Temperature
    plt.subplot(4, 1, 4)
    plt.plot(t[::down2], timeTher['cap']['C1'][::down2])
    plt.ylabel("$\Theta(t)$ (°C)")
    plt.title('Time-domain Thermal DC-Link Capacitor')
    plt.xlabel('time in (sec)')
    plt.grid('on')
    plt.show()

    ###################################################################################################################
    # MSG OUT
    ###################################################################################################################
    print("END: Plotting Stationary Waveforms")
