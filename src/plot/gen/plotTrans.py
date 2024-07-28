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
    plt.plot(t[::down2], timeAc['i_a'][::down2])
    plt.ylabel("$i_{a}(t)$ (A)")
    plt.title('Time-domain Currents AC-Side')
    plt.xlabel('time in (sec)')
    pl.legend(["$i_{a}$"], loc='upper right')
    plt.grid('on')

    # Voltage
    plt.subplot(2, 2, 2)
    plt.plot(t[::down2], timeAc['v_a'][::down2], t[::down2], timeAc['v_a0'][::down2], t[::down2],
             timeAc['v_a_out'][::down2], t[::down2], timeSw['e_a'][::down2])
    plt.ylabel("$v_{a}(t)$ (V)")
    plt.title('Time-domain Voltages AC-Side')
    plt.xlabel('time in (sec)')
    plt.legend(["$v_{a}(t)$", "$v_{a0}(t)$", "$v_{a,out}(t)$", "$e_{a}(t)$"], loc='upper right')
    plt.grid('on')

    # ------------------------------------------
    # DC-Link
    # ------------------------------------------
    # Current
    plt.subplot(2, 2, 3)
    plt.plot(t[::down2], timeDc['i_dc'][::down2], t[::down2],
             np.mean(timeDc['i_dc']) * np.ones(np.size(timeDc['i_dc'][::down2])), '--')
    plt.ylabel("$i_{dc}(t)$ (A)")
    plt.title('Time-domain Currents DC-Side')
    plt.xlabel('time in (sec)')
    plt.legend(["$i_{dc}$", "$I_{dc,avg}$"], loc='upper right')
    plt.grid('on')

    # Voltage
    plt.subplot(2, 2, 4)
    plt.plot(t[::down2], timeDc['v_in'][::down2], t[::down2], timeDc['v_dc'][::down2], t[::down2],
             np.mean(timeDc['v_dc']) * np.ones(np.size(timeDc['v_dc'][::down2])), '--')
    plt.ylabel("$v_{dc}(t)$ (V)")
    plt.title('Time-domain Voltages DC-Side')
    plt.xlabel('time in (sec)')
    plt.legend(["$v_{in}$", "$v_{dc}$", "$V_{dc,avg}$"], loc='upper right')
    plt.grid('on')

    # ==============================================================================
    # Time-domain Transient
    # ==============================================================================
    plt.figure()
    txt = "Time domain switches for PWM control with: " + "$V_{dc}$=" + str(Vdc) + "V, " + "$M_{i}$=" + str(Mi) + ", $\phi_{RL}=$" + str(int(math.degrees(angZ))) + "deg" + ", $\phi_{E}=$" + str(int(phiE)) + "deg" + ", $\phi_{VI}=$" + str(int(math.degrees(phi))) + "deg"
    plt.suptitle(txt, size=18)
    plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)

    # ------------------------------------------
    # Switches Leg 1
    # ------------------------------------------
    plt.subplot(411)
    plt.plot(tel[::down2], timeElec['sw']['S1']['v_T'][::down2], 'b', tel[::down2], timeElec['sw']['S2']['v_T'][::down2], 'r', tel[::down2],
             timeElec['sw']['S1']['v_D'][::down2], 'b--', tel[::down2], timeElec['sw']['S2']['v_D'][::down2], 'r--')
    plt.title('Voltages Transistors and Diodes (A)')
    plt.ylabel('Voltage in (V)')
    plt.xlabel("")
    plt.legend(['T1', 'T2', 'D1', 'D2'])
    plt.grid('on')

    plt.subplot(412)
    plt.plot(tel[::down2], timeElec['sw']['S1']['i_T'][::down2], 'b', tel[::down2], timeElec['sw']['S2']['i_T'][::down2], 'r', tel[::down2],
             timeElec['sw']['S1']['i_D'][::down2], 'b--', tel[::down2], timeElec['sw']['S2']['i_D'][::down2], 'r--')
    plt.title('Currents Transistors and Diodes (A)')
    plt.ylabel('Current in (A)')
    plt.xlabel("")
    plt.legend(['T1', 'T2', 'D1', 'D2'])
    plt.grid('on')

    plt.subplot(413)
    plt.plot(tel[::down2], timeLoss['sw']['S1']['p_T_c'][::down2], 'b', tel[::down2], timeLoss['sw']['S2']['p_T_c'][::down2], 'r', tel[::down2],
             timeLoss['sw']['S1']['p_D_c'][::down2], 'b--', tel[::down2], timeLoss['sw']['S2']['p_D_c'][::down2], 'r--')
    plt.title('Conduction Losses Transistors and Diodes (A)')
    plt.ylabel('Power in (W)')
    plt.xlabel("")
    plt.legend(['T1', 'T2', 'D1', 'D2'])
    plt.grid('on')

    plt.subplot(414)
    plt.plot(tel[::down2], timeLoss['sw']['S1']['p_T_s'][::down2], 'b', tel[::down2], timeLoss['sw']['S2']['p_T_s'][::down2], 'r', tel[::down2],
             timeLoss['sw']['S1']['p_D_s'][::down2], 'b--', tel[::down2], timeLoss['sw']['S2']['p_D_s'][::down2], 'r--')
    plt.title('Switching Losses Transistors and Diodes (A)')
    plt.ylabel('Power in (W)')
    plt.xlabel('time in (sec)')
    plt.legend(['T1', 'T2', 'D1', 'D2'])
    plt.grid('on')

    # ------------------------------------------
    # Capacitor
    # ------------------------------------------
    plt.figure()
    plt.subplot(4, 1, 1)
    plt.plot(tel[::down2], timeElec['cap']['C1']['v_c'][::down2])
    plt.title('Voltage Capacitor')
    plt.ylabel('Voltage in (V)')
    plt.xlabel("")
    plt.grid('on')

    plt.subplot(4, 1, 2)
    plt.plot(tel[::down2], timeElec['cap']['C1']['i_c'][::down2])
    plt.title('Current Capacitor')
    plt.ylabel('Current in (A)')
    plt.xlabel("")
    plt.grid('on')

    plt.subplot(4, 1, 3)
    plt.plot(tel[::down2], timeLoss['cap']['C1']['p_L'][::down2])
    plt.title('Losses Capacitor')
    plt.ylabel('Power in (W)')
    plt.xlabel("")
    plt.grid('on')

    plt.subplot(4, 1, 4)
    plt.plot(tel[::down2], timeTher['cap']['C1'][::down2])
    plt.title('Thermal Capacitor')
    plt.ylabel('Temperature in (°C)')
    plt.grid('on')
    plt.xlabel('time in (sec)')

    # ------------------------------------------
    # Thermal
    # ------------------------------------------
    plt.figure()
    plt.subplot(211)
    plt.plot(tel[::down2], timeLoss['sw']['S1']['p_L'][::down2], 'r', tel[::down2], timeLoss['sw']['S2']['p_L'][::down2], 'b')
    plt.title('Losses Switches (A)')
    plt.ylabel('Power in (W)')
    plt.xlabel("")
    plt.legend(['S1', 'S2'])
    plt.grid('on')

    plt.subplot(212)
    plt.plot(tel[::down2], timeTher['sw']['T1'][::down2], color='b', linestyle='--')
    plt.plot(tel[::down2], timeTher['sw']['D1'][::down2], color='b', linestyle=':')
    plt.plot(tel[::down2], timeTher['sw']['C1'][::down2], color='b', linestyle='-')
    plt.plot(tel[::down2], timeTher['sw']['T2'][::down2], color='r', linestyle='--')
    plt.plot(tel[::down2], timeTher['sw']['D2'][::down2], color='r', linestyle=':')
    plt.plot(tel[::down2], timeTher['sw']['C2'][::down2], color='r', linestyle='-')
    plt.plot(tel[::down2], Ta * np.ones(np.size(tel[::down2])), color='k', linestyle='-')
    plt.title('Thermal Switches (A)')
    plt.ylabel('Temperature in (°C)')
    plt.xlabel('time in (sec)')
    plt.legend(['T1_j', 'T1_d', 'T1_c', 'T2_j', 'T2_d', 'T2_c', 'Ta'])
    plt.grid('on')
    plt.show()

    ###################################################################################################################
    # MSG OUT
    ###################################################################################################################
    print("END: Plotting Transient Waveforms")
