#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         plotSweep_B4
# Date:         14.08.2023
# Author:       Dr. Pascal A. Schirmer
# Version:      V.0.2
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
from src.general.helpFnc import thd

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
def plotSweep_B4(time, freq, sweep, setupPara, setupData, setupTopo, setupExp):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("START: Plotting B4")

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

    # ------------------------------------------
    # Frequency
    # ------------------------------------------
    freqSw = freq['Sw']
    freqAc = freq['Ac']
    freqDc = freq['Dc']

    # ------------------------------------------
    # Sweep
    # ------------------------------------------
    distAc = sweep['Ac']
    distDc = sweep['Dc']

    # ==============================================================================
    # Parameters
    # ==============================================================================
    fel = setupTopo['fel']
    fs = setupPara['PWM']['fs']
    fsim = setupExp['fsim']
    Q = int(fs / fel)
    R = setupTopo['R']
    L = setupTopo['L']
    Mi = setupData['stat']['Mi']
    Vdc = setupData['stat']['Vdc']
    phiE = setupTopo['phiE']
    down = int(setupData['stat']['cyc']) - 2
    down2 = int(fsim/fs/200)
    if down2 < 1:
        down2 = 1

    # ==============================================================================
    # Variables
    # ==============================================================================
    t = timeSw['t'].values
    M_i = sweep['Mi']
    f = fsim * np.linspace(0, 0.5, int(len(t) / 2)) / fel

    ###################################################################################################################
    # Pre-Processing
    ###################################################################################################################
    # ==============================================================================
    # Start and End Plotting
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

    # ==============================================================================
    # Load angle Total
    # ==============================================================================
    Y = fft(timeAc['v_a'])
    phiV = np.angle(Y)[1]
    Y = fft(timeAc['i_a'])
    phiI = np.angle(Y)[1]
    phi = phiV - phiI
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
        setupPara['PWM']['samp']) + ", Update: " + str(setupPara['PWM']['upd']) + " and Edge Trigger: " + str(
        setupPara['PWM']['tri'])
    pl.suptitle(txt, size=18)
    pl.subplots_adjust(hspace=0.35, wspace=0.20, left=0.075, right=0.925, top=0.90, bottom=0.075)

    # ------------------------------------------
    # Modulation
    # ------------------------------------------
    ax1 = pl.subplot(gs[0, :])
    ax2 = ax1.twinx()
    ax1.plot(t[::down2], timeSw['cA'][::down2], 'tab:blue')
    ax2.plot(t[::down2], timeSw['v_a_ref'][::down2] / (Vdc / 2), color='tab:blue')
    ax1.plot(t[::down2], timeSw['cB'][::down2], 'tab:orange')
    ax2.plot(t[::down2], timeSw['v_b_ref'][::down2] / (Vdc / 2), color='tab:orange')
    pl.title('Carrier and Reference Waveform')
    pl.xlabel('t in (sec)')
    ax1.set_ylabel('c(t)', color='black')
    ax2.set_ylabel('v(t) (p.u)', color='black')
    pl.legend(["$c_{a}$", "$v_{a}^{*}$", "$c_{b}$", "$v_{b}^{*}$"], loc='upper right')
    pl.grid('on')

    # ------------------------------------------
    # Switching Waveform
    # ------------------------------------------
    # Time-Domain
    pl.subplot(gs[1, 0])
    pl.plot(t[::down2], timeSw['sA'][::down2])
    pl.plot(t[::down2], timeSw['sB'][::down2])
    pl.title('Time-domain Switching Functions')
    pl.xlabel('t in (sec)')
    pl.ylabel("$s_{x}(t)$ (p.u)")
    pl.legend(["$s_{a}$", "$s_{b}$"], loc='upper right')
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
    pl.plot(t[::down2], timeSw['xBs'][::down2])
    pl.title('Time-domain Sampled References')
    pl.xlabel('t in (sec)')
    pl.ylabel("$x_{x}(t)$ (p.u)")
    pl.legend(["$x_{a}$", "$x_{b}$"], loc='upper right')
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
    txt = "Currents B4 Bridge for PWM Control with: " + "$V_{dc}$=" + str(Vdc) + "V, " + "$M_{i}$=" + str(
        Mi) + "$ ,Q$=" + str(Q) + ", $\phi_{RL}=$" + str(int(math.degrees(angZ))) + "deg" + ", $\phi_{E}=$" + str(
        int(phiE)) + "deg" + ", $\phi_{VI}=$" + str(int(math.degrees(phi))) + "deg"
    plt.suptitle(txt, size=18)
    plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)

    # ------------------------------------------
    # Phase
    # ------------------------------------------
    # Time
    plt.subplot(2, 3, 1)
    plt.plot(t[::down2], timeAc['i_a'][::down2])
    plt.ylabel("$i_{a}(t)$ (A)")
    plt.title('Time-domain Currents AC-Side')
    plt.xlabel('time in (sec)')
    plt.grid('on')

    # Frequency
    ax = plt.subplot(2, 3, 2)
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

    # Modulation
    plt.subplot(233)
    plt.plot(M_i, distAc['num']['I_a_thd'])
    if setupExp['plot'] >= 2:
        plt.plot(M_i, distAc['ana']['I_a_thd'], 'tab:blue', linestyle="", marker="o")
        plt.legend(['Numerical', 'Analytical'])
    plt.ylim(0, )
    plt.title('Distortion Current AC-Side')
    plt.ylabel("$I_{a,rms}^{THD}$ (A)")
    plt.xlabel("$M_{i}$ in (p.u)")
    plt.grid('on')

    # ------------------------------------------
    # DC-Link
    # ------------------------------------------
    # Time
    plt.subplot(2, 3, 4)
    plt.plot(t[::down2], timeDc['i_dc'][::down2], t[::down2], np.mean(timeDc['i_dc']) * np.ones(np.size(timeDc['i_dc'][::down2])), '--')
    plt.ylabel("$i_{dc}(t)$ (A)")
    plt.title('Time-domain Currents DC-Side')
    plt.xlabel('time in (sec)')
    plt.legend(["$i_{dc}$", "$I_{dc,avg}$"], loc='upper right')
    plt.grid('on')

    # Frequency
    ax = plt.subplot(2, 3, 5)
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

    # Modulation
    plt.subplot(236)
    plt.plot(M_i, distDc['num']['I_dc_thd'])
    # if setupExp['plot'] >= 2:
    #    plt.plot(M_i, distDc['ana']['I_dc_thd'], 'tab:blue', linestyle="", marker="o")
    #    plt.legend(['Numerical', 'Analytical'])
    plt.ylim(0, )
    plt.title('Distortion Current DC-Side')
    plt.ylabel("$I_{dc,rms}^{THD}$ (A)")
    plt.xlabel("$M_{i}$ in (p.u)")
    plt.grid('on')

    # ==============================================================================
    # Voltage
    # ==============================================================================
    plt.figure()
    txt = "Voltages B4 Bridge for PWM Control with: " + "$V_{dc}$=" + str(Vdc) + "V, " + "$M_{i}$=" + str(
        Mi) + ", $\phi_{RL}=$" + str(int(math.degrees(angZ))) + "deg" + ", $\phi_{E}=$" + str(
        int(phiE)) + "deg" + ", $\phi_{VI}=$" + str(int(math.degrees(phi))) + "deg"
    plt.suptitle(txt, size=18)
    plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)

    # ------------------------------------------
    # Phase
    # ------------------------------------------
    # Time
    plt.subplot(2, 3, 1)
    plt.plot(t[::down2], timeAc['v_a'][::down2], t[::down2], timeAc['v_a0'][::down2], t[::down2], timeSw['e'][::down2])
    plt.ylabel("$v_{ab}(t)$ (V)")
    plt.title('Time-domain Voltages AC-Side')
    plt.xlabel('time in (sec)')
    plt.legend(["$v_{ab}(t)$", "$v_{a0}(t)$", "$e(t)$"], loc='upper right')
    plt.grid('on')

    # Frequency
    ax = plt.subplot(2, 3, 2)
    plt.stem(f[::down][0:50], freqAc['V_a'][::down][0:50])
    plt.ylabel("$V_{ab}(f)$ (V)")
    plt.xlim(0, 50)
    plt.title('Frequency-domain Voltages AC-Side')
    plt.xlabel("$f/f_{1}$ (Hz/Hz)")
    plt.yscale('log')
    plt.ylim(0.1 / OoM(max(freqAc['V_a'])), )
    txt = "THD=" + str(round(V_ph_thd * 100, 2)) + "%"
    plt.text(0.75, 0.90, txt, bbox=dict(facecolor='tab:blue', alpha=0.5), transform=ax.transAxes)
    plt.grid('on')

    # Modulation
    plt.subplot(233)
    plt.plot(M_i, distAc['num']['V_a_thd'])
    if setupExp['plot'] >= 2:
        plt.plot(M_i, distAc['ana']['V_a_thd'], 'tab:blue', linestyle="", marker="o")
        plt.legend(['Numerical', 'Analytical'])
    plt.title('Distortion Voltage AC-Side')
    plt.ylabel("$V_{ab,rms}^{THD}$ (V)")
    plt.xlabel("$M_{i}$ in (p.u)")
    plt.grid('on')

    # ------------------------------------------
    # DC-Link
    # ------------------------------------------
    # Time
    plt.subplot(2, 3, 4)
    plt.plot(t[::down2], timeDc['v_in'][::down2], t[::down2], timeDc['v_dc'][::down2], t[::down2], np.mean(timeDc['v_dc']) * np.ones(np.size(timeDc['v_dc'][::down2])), '--')
    plt.ylabel("$v_{dc}(t)$ (V)")
    plt.title('Time-domain Voltages DC-Side')
    plt.xlabel('time in (sec)')
    plt.legend(["$v_{in}$", "$v_{dc}$", "$V_{dc,avg}$"], loc='upper right')
    plt.grid('on')

    # Frequency
    ax = plt.subplot(2, 3, 5)
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

    # Modulation
    plt.subplot(236)
    plt.plot(M_i, distDc['num']['V_dc_thd'])
    # if setupExp['plot'] >= 2:
    #    plt.plot(M_i, distDc['ana']['V_dc_thd'], 'tab:blue', linestyle="", marker="o")
    #    plt.legend(['Numerical', 'Analytical'])
    plt.title('Distortion Voltage DC-Side')
    plt.ylabel("$V_{dc,rms}^{THD}$ (V)")
    plt.xlabel("$M_{i}$ in (p.u)")
    plt.grid('on')
    plt.show()

    ###################################################################################################################
    # Post-Processing
    ###################################################################################################################

    ###################################################################################################################
    # MSG OUT
    ###################################################################################################################
    print("END: Plotting B4")

    ###################################################################################################################
    # Return
    ###################################################################################################################
