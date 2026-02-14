#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         plotTrans_DAB
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
from src.general.helpFnc import OoM

# ==============================================================================
# External
# ==============================================================================
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import matplotlib.gridspec as gridspec
import matplotlib
matplotlib.use('TkAgg')


#######################################################################################################################
# Function
#######################################################################################################################
def plotTrans_DAB(time, freq, setup):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("START: Plotting Transient DAB Waveforms")

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
    phiDAB = setup['Dat']['stat']['PhiDAB']
    Ta = setup['Dat']['trans']['Tc']
    down2 = int(fsim / fs / 200)
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
    """
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
    """

    # ==============================================================================
    # Current/Voltage
    # ==============================================================================
    plt.figure()
    txt = "Switching Functions DAB for PWM Control with: " \
          + "$V_{dc}$=" + str(Vdc) + "V, " \
          + "$M_{i}$=" + str(Mi) + "$ ,Q$=" + str(Q) \
          + ", Phi=" + f"{math.degrees(phiDAB):.2f}"
    plt.suptitle(txt, size=18)
    plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)

    # ------------------------------------------
    # Phase
    # ------------------------------------------
    # Current
    plt.subplot(2, 2, 1)
    plt.plot(t_plot[::down2], timeAc['i_ac_pri'][::down2])
    plt.plot(t_plot[::down2], timeAc['i_ac_sec'][::down2])
    plt.ylabel("$i_{ac}(t)$ (A)")
    plt.title('Time-domain Currents AC')
    plt.xlabel(t_label)
    pl.legend(["$i_{ac,pri}$", "$i_{ac,sec}$"], loc='upper right')
    plt.grid(True)

    # Voltage
    plt.subplot(2, 2, 2)
    plt.plot(t_plot[::down2], timeAc['v_ac_pri'][::down2])
    plt.plot(t_plot[::down2], timeAc['v_ac_sec'][::down2])
    plt.ylabel("$v_{ac}(t)$ (V)")
    plt.title('Time-domain Voltages AC')
    plt.xlabel(t_label)
    pl.legend(["$v_{ac,pri}$", "$v_{ac,sec}$"], loc='upper right')
    plt.grid(True)

    # ------------------------------------------
    # DC-Link
    # ------------------------------------------
    # Current
    plt.subplot(2, 2, 3)
    plt.plot(t_plot[::down2], timeDc['i_dc_pri'][::down2])
    plt.plot(t_plot[::down2], timeDc['i_dc_sec'][::down2])
    plt.ylabel("$i_{dc}(t)$ (A)")
    plt.title('Time-domain Currents DC')
    plt.xlabel(t_label)
    pl.legend(["$i_{dc,pri}$", "$i_{dc,sec}$"], loc='upper right')
    plt.grid(True)

    # Voltage
    plt.subplot(2, 2, 4)
    plt.plot(t_plot[::down2], timeDc['v_dc_pri'][::down2])
    plt.plot(t_plot[::down2], timeDc['v_dc_sec'][::down2])
    plt.ylabel("$v_{dc}(t)$ (V)")
    plt.title('Time-domain Voltages DC')
    plt.xlabel(t_label)
    pl.legend(["$v_{dc,pri}$", "$v_{dc,sec}$"], loc='upper right')
    plt.grid(True)

    # ==============================================================================
    # Time-domain Transient
    # ==============================================================================
    fig, axs = plt.subplots(4, 1, sharex=True)
    txt = "Switching Functions DAB for PWM Control with: " \
          + "$V_{dc}$=" + str(Vdc) + "V, " \
          + "$M_{i}$=" + str(Mi) + "$ ,Q$=" + str(Q) \
          + ", Phi=" + f"{math.degrees(phiDAB):.2f}"
    plt.suptitle(txt, size=18)
    plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)

    # ------------------------------------------
    # Switches Leg 1
    # ------------------------------------------
    # Switches Leg 1 - Voltages
    axs[0].plot(tel[::down2], timeElec['sw']['S1']['v_T'][::down2], 'b',
                tel[::down2], timeElec['sw']['S5']['v_T'][::down2], 'r',
                tel[::down2], timeElec['sw']['S1']['v_D'][::down2], 'b--',
                tel[::down2], timeElec['sw']['S5']['v_D'][::down2], 'r--')
    axs[0].set_title('Voltages Transistors and Diodes (Pri/Sec)')
    axs[0].set_ylabel('Voltage (V)')
    axs[0].legend(['T1', 'T5', 'D1', 'D5'])
    axs[0].grid(True)

    # Switches Leg 1 - Currents
    axs[1].plot(tel[::down2], timeElec['sw']['S1']['i_T'][::down2], 'b',
                tel[::down2], timeElec['sw']['S5']['i_T'][::down2], 'r',
                tel[::down2], timeElec['sw']['S1']['i_D'][::down2], 'b--',
                tel[::down2], timeElec['sw']['S5']['i_D'][::down2], 'r--')
    axs[1].set_title('Currents Transistors and Diodes (Pri/Sec)')
    axs[1].set_ylabel('Current (A)')
    axs[1].legend(['T1', 'T5', 'D1', 'D5'])
    axs[1].grid(True)

    # Switches Leg 1 - Conduction Losses
    axs[2].plot(tel[::down2], timeLoss['sw']['S1']['p_T_c'][::down2], 'b',
                tel[::down2], timeLoss['sw']['S5']['p_T_c'][::down2], 'r',
                tel[::down2], timeLoss['sw']['S1']['p_D_c'][::down2], 'b--',
                tel[::down2], timeLoss['sw']['S5']['p_D_c'][::down2], 'r--')
    axs[2].set_title('Conduction Losses Transistors and Diodes (Pri/Sec)')
    axs[2].set_ylabel('Power (W)')
    axs[2].legend(['T1', 'T5', 'D1', 'D5'])
    axs[2].grid(True)

    # Switches Leg 1 - Switching Losses
    axs[3].plot(tel[::down2], timeLoss['sw']['S1']['p_T_s'][::down2], 'b',
                tel[::down2], timeLoss['sw']['S5']['p_T_s'][::down2], 'r',
                tel[::down2], timeLoss['sw']['S1']['p_D_s'][::down2], 'b--',
                tel[::down2], timeLoss['sw']['S5']['p_D_s'][::down2], 'r--')
    axs[3].set_title('Switching Losses Transistors and Diodes (Pri/Sec)')
    axs[3].set_ylabel('Power (W)')
    axs[3].set_xlabel('Time (sec)')
    axs[3].legend(['T1', 'T5', 'D1', 'D5'])
    axs[3].grid(True)
    plt.show()

    # ------------------------------------------
    # Capacitor
    # ------------------------------------------
    fig, axs = plt.subplots(4, 2, sharex=True)
    txt = "Switching Functions DAB for PWM Control with: " \
          + "$V_{dc}$=" + str(Vdc) + "V, " \
          + "$M_{i}$=" + str(Mi) + "$ ,Q$=" + str(Q) \
          + ", Phi=" + f"{math.degrees(phiDAB):.2f}"
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
    txt = "Switching Functions DAB for PWM Control with: " \
          + "$V_{dc}$=" + str(Vdc) + "V, " \
          + "$M_{i}$=" + str(Mi) + "$ ,Q$=" + str(Q) \
          + ", Phi=" + f"{math.degrees(phiDAB):.2f}"
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
    print("END: Plotting Transient DAB Waveforms")
