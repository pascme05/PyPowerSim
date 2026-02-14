#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         plotTrans
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
def plotTrans(time, freq, setup):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("START: Plotting Transient Waveforms")

    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
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
    down = int(setup['Dat']['stat']['cyc']) - 1
    Ta = setup['Dat']['trans']['Tc']
    down2 = int(fsim/fs/200)
    if down2 < 1:
        down2 = 1

    # ==============================================================================
    # Variables
    # ==============================================================================
    t = timeSw['t'].values
    tel = time['t'] * int(len(time['t']) / len(t))
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
    start = 0
    ende = int((len(t) - 1) / K) + 1

    # ------------------------------------------
    # Change time
    # ------------------------------------------
    t = t[start:ende]
    timeSw = timeSw[:][start:ende]
    for c1 in timeAc:
        timeAc[c1] = timeAc[c1][start:ende]
    for c1 in timeDc:
        timeDc[c1] = timeDc[c1][start:ende]

    # ------------------------------------------
    # Plot Axis
    # ------------------------------------------
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
    ax1.plot(t_plot[::down2], timeSw['c'][::down2], 'black')
    ax2.plot(t_plot[::down2], timeSw['v_a_ref'][::down2] / (Vdc / 2), color='tab:blue')
    ax2.plot(t_plot[::down2], timeSw['xA'][::down2], '--', color='tab:blue')
    ax2.plot(t_plot[::down2], timeSw['n0'][::down2], '-.', color='tab:blue')
    pl.title('Carrier and Reference Waveform')
    pl.xlabel(t_label)
    ax1.set_ylabel('c(t)/s(t)', color='black')
    ax2.set_ylabel('v(t) (p.u)', color='tab:blue')
    pl.legend(["$c_{a}$", "$v_{a}^{*}$", "$v_{a0}^{*}$", "$v_{n0}^{*}$"], loc='upper right')
    pl.grid(True)

    # ------------------------------------------
    # Switching Waveform
    # ------------------------------------------
    # Time-Domain
    pl.subplot(gs[1, 0])
    pl.plot(t_plot[::down2], timeSw['sA'][::down2])
    pl.title('Time-domain Switching Functions')
    pl.xlabel('t in (sec)')
    pl.ylabel("$s_{a}(t)$ (p.u)")
    pl.legend(["$s_{a}$"], loc='upper right')
    pl.grid(True)

    # Freq-Domain
    pl.subplot(gs[2, 0])
    pl.stem(f[::down][0:50], freqSw['Sa'][::down][0:50])
    pl.xlim(0, 50)
    pl.title('Frequency-domain Switching Function')
    pl.xlabel("$f/f_{1}$ (Hz/Hz)")
    pl.ylabel("$S_{a}(f)$ (p.u)")
    pl.yscale('log')
    pl.ylim(0.0001, )
    pl.grid(True)

    # ------------------------------------------
    # Reference Waveform
    # ------------------------------------------
    # Time-Domain
    pl.subplot(gs[1, 1])
    pl.plot(t_plot[::down2], timeSw['xAs'][::down2])
    pl.title('Time-domain Sampled References')
    pl.xlabel('t in (sec)')
    pl.ylabel("$x_{a}(t)$ (p.u)")
    pl.legend(["$x_{a}$"], loc='upper right')
    pl.grid(True)

    # Freq-Domain
    pl.subplot(gs[2, 1])
    pl.stem(f[::down][0:50], freqSw['Xas'][::down][0:50])
    pl.xlim(0, 50)
    pl.title('Frequency-domain Sampled Reference')
    pl.xlabel("$f/f_{1}$ (Hz/Hz)")
    pl.ylabel("$X_{a}(f)$ (p.u)")
    pl.yscale('log')
    pl.ylim(0.0001, )
    pl.grid(True)

    # ==============================================================================
    # Current/Voltage
    # ==============================================================================
    plt.figure()
    txt = "Currents and Voltages for PWM control with: " + "$V_{dc}$=" + str(Vdc) + "V, " + "$M_{i}$=" + str(Mi) + "$ ,Q$=" + str(Q) + ", $\phi_{RL}=$" + str(int(math.degrees(angZ))) + "deg" + ", $\phi_{E}=$" + str(int(phiE)) + "deg" + ", $\phi_{VI}=$" + str(int(math.degrees(phi))) + "deg"
    plt.suptitle(txt, size=18)
    plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)

    # ------------------------------------------
    # Phase
    # ------------------------------------------
    # Current
    plt.subplot(2, 2, 1)
    plt.plot(t_plot[::down2], timeAc['i_a'][::down2])
    plt.ylabel("$i_{a}(t)$ (A)")
    plt.title('Time-domain Currents AC-Side')
    plt.xlabel(t_label)
    pl.legend(["$i_{a}$"], loc='upper right')
    plt.grid(True)

    # Voltage
    plt.subplot(2, 2, 2)
    plt.plot(t_plot[::down2], timeAc['v_a'][::down2], t_plot[::down2], timeAc['v_a0'][::down2], t_plot[::down2],
             timeAc['v_a_out'][::down2], t_plot[::down2], timeSw['e_a'][::down2])
    plt.ylabel("$v_{a}(t)$ (V)")
    plt.title('Time-domain Voltages AC-Side')
    plt.xlabel(t_label)
    plt.legend(["$v_{a}(t)$", "$v_{a0}(t)$", "$v_{a,out}(t)$", "$e_{a}(t)$"], loc='upper right')
    plt.grid(True)

    # ------------------------------------------
    # DC-Link
    # ------------------------------------------
    # Current
    plt.subplot(2, 2, 3)
    plt.plot(t_plot[::down2], timeDc['i_dc'][::down2], t_plot[::down2],
             np.mean(timeDc['i_dc']) * np.ones(np.size(timeDc['i_dc'][::down2])), '--')
    plt.ylabel("$i_{dc}(t)$ (A)")
    plt.title('Time-domain Currents DC-Side')
    plt.xlabel(t_label)
    plt.legend(["$i_{dc}$", "$I_{dc,avg}$"], loc='upper right')
    plt.grid(True)

    # Voltage
    plt.subplot(2, 2, 4)
    plt.plot(t_plot[::down2], timeDc['v_in'][::down2], t_plot[::down2], timeDc['v_dc'][::down2], t_plot[::down2],
             np.mean(timeDc['v_dc']) * np.ones(np.size(timeDc['v_dc'][::down2])), '--')
    plt.ylabel("$v_{dc}(t)$ (V)")
    plt.title('Time-domain Voltages DC-Side')
    plt.xlabel(t_label)
    plt.legend(["$v_{in}$", "$v_{dc}$", "$V_{dc,avg}$"], loc='upper right')
    plt.grid(True)

    # ==============================================================================
    # Time-domain Transient
    # ==============================================================================
    fig, axs = plt.subplots(4, 1, sharex=True)
    txt = ("Time domain switches for PWM control with: $V_{dc}$=" + str(Vdc) + "V, $M_{i}$=" + str(Mi) +
           ", $\phi_{RL}=$" + str(int(math.degrees(angZ))) + "deg, $\phi_{E}=$" + str(int(phiE)) +
           "deg, $\phi_{VI}=$" + str(int(math.degrees(phi))) + "deg")
    plt.suptitle(txt, size=18)
    plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)

    # ------------------------------------------
    # Switches Leg 1
    # ------------------------------------------
    # Switches Leg 1 - Voltages
    axs[0].plot(tel[::down2], timeElec['sw']['S1']['v_T'][::down2], 'b',
                tel[::down2], timeElec['sw']['S2']['v_T'][::down2], 'r',
                tel[::down2], timeElec['sw']['S1']['v_D'][::down2], 'b--',
                tel[::down2], timeElec['sw']['S2']['v_D'][::down2], 'r--')
    axs[0].set_title('Voltages Transistors and Diodes (A)')
    axs[0].set_ylabel('Voltage (V)')
    axs[0].legend(['T1', 'T2', 'D1', 'D2'])
    axs[0].grid(True)

    # Switches Leg 1 - Currents
    axs[1].plot(tel[::down2], timeElec['sw']['S1']['i_T'][::down2], 'b',
                tel[::down2], timeElec['sw']['S2']['i_T'][::down2], 'r',
                tel[::down2], timeElec['sw']['S1']['i_D'][::down2], 'b--',
                tel[::down2], timeElec['sw']['S2']['i_D'][::down2], 'r--')
    axs[1].set_title('Currents Transistors and Diodes (A)')
    axs[1].set_ylabel('Current (A)')
    axs[1].legend(['T1', 'T2', 'D1', 'D2'])
    axs[1].grid(True)

    # Switches Leg 1 - Conduction Losses
    axs[2].plot(tel[::down2], timeLoss['sw']['S1']['p_T_c'][::down2], 'b',
                tel[::down2], timeLoss['sw']['S2']['p_T_c'][::down2], 'r',
                tel[::down2], timeLoss['sw']['S1']['p_D_c'][::down2], 'b--',
                tel[::down2], timeLoss['sw']['S2']['p_D_c'][::down2], 'r--')
    axs[2].set_title('Conduction Losses Transistors and Diodes (A)')
    axs[2].set_ylabel('Power (W)')
    axs[2].legend(['T1', 'T2', 'D1', 'D2'])
    axs[2].grid(True)

    # Switches Leg 1 - Switching Losses
    axs[3].plot(tel[::down2], timeLoss['sw']['S1']['p_T_s'][::down2], 'b',
                tel[::down2], timeLoss['sw']['S2']['p_T_s'][::down2], 'r',
                tel[::down2], timeLoss['sw']['S1']['p_D_s'][::down2], 'b--',
                tel[::down2], timeLoss['sw']['S2']['p_D_s'][::down2], 'r--')
    axs[3].set_title('Switching Losses Transistors and Diodes (A)')
    axs[3].set_ylabel('Power (W)')
    axs[3].set_xlabel('Time (sec)')
    axs[3].legend(['T1', 'T2', 'D1', 'D2'])
    axs[3].grid(True)

    # ------------------------------------------
    # Capacitor
    # ------------------------------------------
    fig, axs = plt.subplots(4, 1, sharex=True)
    txt = ("Time domain capacitor for PWM control with: $V_{dc}$=" + str(Vdc) + "V, $M_{i}$=" + str(Mi) +
           ", $\phi_{RL}=$" + str(int(math.degrees(angZ))) + "deg, $\phi_{E}=$" + str(int(phiE)) +
           "deg, $\phi_{VI}=$" + str(int(math.degrees(phi))) + "deg")
    plt.suptitle(txt, size=18)
    plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)

    # Voltage Capacitor
    axs[0].plot(tel[::down2], timeElec['cap']['C1']['v_c'][::down2])
    axs[0].set_title('Voltage Capacitor')
    axs[0].set_ylabel('Voltage (V)')
    axs[0].grid(True)

    # Current Capacitor
    axs[1].plot(tel[::down2], timeElec['cap']['C1']['i_c'][::down2])
    axs[1].set_title('Current Capacitor')
    axs[1].set_ylabel('Current (A)')
    axs[1].grid(True)

    # Losses Capacitor
    axs[2].plot(tel[::down2], timeLoss['cap']['C1']['p_L'][::down2])
    axs[2].set_title('Losses Capacitor')
    axs[2].set_ylabel('Power (W)')
    axs[2].grid(True)

    # Thermal Capacitor
    axs[3].plot(tel[::down2], timeTher['cap']['C1'][::down2])
    axs[3].set_title('Thermal Capacitor')
    axs[3].set_ylabel('Temperature (°C)')
    axs[3].set_xlabel('Time (sec)')
    axs[3].grid(True)

    # ------------------------------------------
    # Thermal
    # ------------------------------------------
    fig, axs = plt.subplots(2, 1, sharex=True)
    txt = ("Time domain thermal for PWM control with: $V_{dc}$=" + str(Vdc) + "V, $M_{i}$=" + str(Mi) +
           ", $\phi_{RL}=$" + str(int(math.degrees(angZ))) + "deg, $\phi_{E}=$" + str(int(phiE)) +
           "deg, $\phi_{VI}=$" + str(int(math.degrees(phi))) + "deg")
    plt.suptitle(txt, size=18)
    plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)

    # Losses Switches
    axs[0].plot(tel[::down2], timeLoss['sw']['S1']['p_L'][::down2], 'r',
                tel[::down2], timeLoss['sw']['S2']['p_L'][::down2], 'b')
    axs[0].set_title('Losses Switches (A)')
    axs[0].set_ylabel('Power (W)')
    axs[0].legend(['S1', 'S2'])
    axs[0].grid(True)

    # Thermal Switches
    axs[1].plot(tel[::down2], timeTher['sw']['T1'][::down2], color='b', linestyle='--')
    axs[1].plot(tel[::down2], timeTher['sw']['D1'][::down2], color='b', linestyle=':')
    axs[1].plot(tel[::down2], timeTher['sw']['C1'][::down2], color='b', linestyle='-')
    axs[1].plot(tel[::down2], timeTher['sw']['T2'][::down2], color='r', linestyle='--')
    axs[1].plot(tel[::down2], timeTher['sw']['D2'][::down2], color='r', linestyle=':')
    axs[1].plot(tel[::down2], timeTher['sw']['C2'][::down2], color='r', linestyle='-')
    axs[1].plot(tel[::down2], Ta * np.ones(np.size(tel[::down2])), color='k', linestyle='-')
    axs[1].set_title('Thermal Switches (A)')
    axs[1].set_ylabel('Temperature (°C)')
    axs[1].set_xlabel('Time (sec)')
    axs[1].legend(['T1_j', 'T1_d', 'T1_c', 'T2_j', 'T2_d', 'T2_c', 'Ta'])
    axs[1].grid(True)
    plt.show()

    ###################################################################################################################
    # MSG OUT
    ###################################################################################################################
    print("END: Plotting Transient Waveforms")
