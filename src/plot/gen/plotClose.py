#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         plotClose
# Date:         07.08.2024
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
from scipy.fft import fft
import matplotlib
matplotlib.use('TkAgg')


#######################################################################################################################
# Function
#######################################################################################################################
def plotClose(time, freq, setup):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("START: Plotting Closed Loop Waveforms")

    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Init
    # ==============================================================================
    timeSw = time['Sw']
    timeAc = time['Ac']
    timeElec = time['Elec']
    timeLoss = time['Loss']

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
    down2 = int(fsim/fs/200)
    if down2 < 1:
        down2 = 1

    # ==============================================================================
    # Variables
    # ==============================================================================
    tel = time['t'] * setup['Dat']['trans']['tmax']*fel

    ###################################################################################################################
    # Pre-Processing
    ###################################################################################################################
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
    # Phase
    # ==============================================================================
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
    txt = ("Currents and Voltages for PWM control with: $V_{dc}$=" + str(Vdc) + "V, $M_{i}$=" + str(Mi) +
           "$ ,Q$=" + str(Q) + ", $\phi_{RL}=$" + str(int(math.degrees(angZ))) + "deg, $\phi_{E}=$" +
           str(int(phiE)) + "deg, $\phi_{VI}=$" + str(int(math.degrees(phi))) + "deg")
    fig.suptitle(txt, size=18)
    fig.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)

    # Current
    ax1.plot(tel[::down2], timeAc['i_a'][::down2])
    ax1.plot(tel[::down2], timeAc['i_ref']['A'][::down2])
    ax1.set_ylabel("$i_{a}(t)$ (A)")
    ax1.set_title('Time-domain Currents AC-Side')
    ax1.legend(["$i_{a}^{act}$", "$i_{a}^{ref}$"], loc='upper right')
    ax1.grid(True)

    # Current Error
    err = timeAc['i_a'][::down2] - timeAc['i_ref']['A'][::down2]
    ax2.plot(tel[::down2], err)
    ax2.set_ylabel("$i_{err}(t)$ (A)")
    ax2.set_title('Time-domain Error Current AC-Side')
    ax2.legend(["$i_{a}^{err}$"], loc='upper right')
    ax2.grid(True)

    # Voltage
    ax3.plot(tel[::down2], timeAc['v_a'][::down2], label="$v_{a}(t)$")
    ax3.plot(tel[::down2], timeAc['v_a0'][::down2], label="$v_{a0}(t)$")
    ax3.plot(tel[::down2], timeAc['v_a_out'][::down2], label="$v_{a,out}(t)$")
    ax3.set_ylabel("$v_{a}(t)$ (V)")
    ax3.set_title('Time-domain Voltages AC-Side')
    ax3.set_xlabel('Time (sec)')
    ax3.legend(loc='upper right')
    ax3.grid(True)

    # ==============================================================================
    # Time-domain Transient
    # ==============================================================================
    fig, axs = plt.subplots(4, 1, sharex=True)
    txt = "Time domain switches for PWM control with: " + "$V_{dc}$=" + str(Vdc) + "V, " + "$M_{i}$=" + str(
        Mi) + ", $\phi_{RL}=$" + str(int(math.degrees(angZ))) + "deg" + ", $\phi_{E}=$" + str(
        int(phiE)) + "deg" + ", $\phi_{VI}=$" + str(int(math.degrees(phi))) + "deg"
    plt.suptitle(txt, size=18)
    plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)

    # ------------------------------------------
    # Switches Leg 1
    # ------------------------------------------
    # Voltages Transistors and Diodes (A)
    axs[0].plot(tel[::down2], timeElec['sw']['S1']['v_T'][::down2], 'b', tel[::down2],
                timeElec['sw']['S2']['v_T'][::down2], 'r',
                tel[::down2], timeElec['sw']['S1']['v_D'][::down2], 'b--', tel[::down2],
                timeElec['sw']['S2']['v_D'][::down2], 'r--')
    axs[0].set_title('Voltages Transistors and Diodes (A)')
    axs[0].set_ylabel('Voltage (V)')
    axs[0].legend(['T1', 'T2', 'D1', 'D2'])
    axs[0].grid(True)

    # Currents Transistors and Diodes (A)
    axs[1].plot(tel[::down2], timeElec['sw']['S1']['i_T'][::down2], 'b', tel[::down2],
                timeElec['sw']['S2']['i_T'][::down2], 'r',
                tel[::down2], timeElec['sw']['S1']['i_D'][::down2], 'b--', tel[::down2],
                timeElec['sw']['S2']['i_D'][::down2], 'r--')
    axs[1].set_title('Currents Transistors and Diodes (A)')
    axs[1].set_ylabel('Current (A)')
    axs[1].legend(['T1', 'T2', 'D1', 'D2'])
    axs[1].grid(True)

    # Conduction Losses Transistors and Diodes (A)
    axs[2].plot(tel[::down2], timeLoss['sw']['S1']['p_T_c'][::down2], 'b', tel[::down2],
                timeLoss['sw']['S2']['p_T_c'][::down2], 'r',
                tel[::down2], timeLoss['sw']['S1']['p_D_c'][::down2], 'b--', tel[::down2],
                timeLoss['sw']['S2']['p_D_c'][::down2], 'r--')
    axs[2].set_title('Conduction Losses Transistors and Diodes (A)')
    axs[2].set_ylabel('Power (W)')
    axs[2].legend(['T1', 'T2', 'D1', 'D2'])
    axs[2].grid(True)

    # Switching Losses Transistors and Diodes (A)
    axs[3].plot(tel[::down2], timeLoss['sw']['S1']['p_T_s'][::down2], 'b', tel[::down2],
                timeLoss['sw']['S2']['p_T_s'][::down2], 'r',
                tel[::down2], timeLoss['sw']['S1']['p_D_s'][::down2], 'b--', tel[::down2],
                timeLoss['sw']['S2']['p_D_s'][::down2], 'r--')
    axs[3].set_title('Switching Losses Transistors and Diodes (A)')
    axs[3].set_ylabel('Power (W)')
    axs[3].set_xlabel('Time (sec)')
    axs[3].legend(['T1', 'T2', 'D1', 'D2'])
    axs[3].grid(True)

    # ------------------------------------------
    # Capacitor
    # ------------------------------------------
    fig, axs = plt.subplots(3, 1, sharex=True)
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
    axs[2].set_xlabel('Time (sec)')
    axs[2].grid(True)
    plt.show()

    ###################################################################################################################
    # MSG OUT
    ###################################################################################################################
    print("END: Plotting Closed Loop Waveforms")
