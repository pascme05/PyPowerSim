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
def plotSweep_DAB(time, freq, sweep, setup):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("START: Plotting DAB Sweep Results")

    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    timeSw = time['Sw']
    timeAc = time['Ac']
    timeDc = time['Dc']
    freqSw = freq['Sw']
    freqAc = freq['Ac']
    freqDc = freq['Dc']
    distAc = sweep['Ac']
    distDc = sweep['Dc']

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
    down = int(setup['Dat']['stat']['cyc']) - 2
    down2 = int(fsim / fs / 200)
    if down2 < 1:
        down2 = 1

    t = timeSw['t'].values
    M_i = sweep['Mi']
    f = fsim * np.linspace(0, 0.5, int(len(t) / 2)) / fel

    ###################################################################################################################
    # Pre-Processing
    ###################################################################################################################
    K = int(np.round((t[-1] - t[0]) * fel))
    start = int((len(t) - 1) / K) * (K - 1)
    ende = len(t)

    t = t[start:ende]
    timeSw = timeSw[:][start:ende]
    for c1 in timeAc:
        timeAc[c1] = timeAc[c1][start:ende]
    for c1 in timeDc:
        timeDc[c1] = timeDc[c1][start:ende]

    ###################################################################################################################
    # Load angle and THD
    ###################################################################################################################
    Y = fft(timeAc['v_a'])
    phiV = np.angle(Y)[1]
    Y = fft(timeAc['i_a'])
    phiI = np.angle(Y)[1]
    phi = phiV - phiI
    while phi > 2 * np.pi:
        phi = phi - 2 * np.pi

    angZ = math.atan2(2 * np.pi * fel * L, R)
    I_ph_thd = thd(timeAc['i_a'], t, (2 / np.sqrt(2)), 1)
    V_ph_thd = thd(timeAc['v_a'], t, (2 / np.sqrt(2)), 1)
    I_dc_thd = thd(timeDc['i_dc'], t, 1, 0)
    V_dc_thd = thd(timeDc['v_dc'], t, 1, 0)

    ###################################################################################################################
    # Switching Functions
    ###################################################################################################################
    gs = gridspec.GridSpec(3, 2)
    pl.figure()
    txt = "DAB Modulation: " + "$M_{i}$=" + str(Mi) + "$ ,Q$=" + str(Q) + ", n=" + str(n_tr)
    pl.suptitle(txt, size=18)
    pl.subplots_adjust(hspace=0.35, wspace=0.20, left=0.075, right=0.925, top=0.90, bottom=0.075)

    ax1 = pl.subplot(gs[0, :])
    ax2 = ax1.twinx()
    ax1.plot(t[::down2], timeSw['c'][::down2], 'black')
    ax2.plot(t[::down2], timeSw['sA'][::down2], color='tab:blue')
    ax2.plot(t[::down2], timeSw['sC'][::down2], color='tab:orange')
    pl.title('Primary and Secondary Switching States')
    pl.xlabel('t in (sec)')
    ax1.set_ylabel('c(t)', color='black')
    ax2.set_ylabel('s(t) (p.u)', color='tab:blue')
    pl.legend(["$c$", "$s_{p}$", "$s_{s}$"], loc='upper right')
    pl.grid('on')

    pl.subplot(gs[1, 0])
    pl.plot(t[::down2], timeSw['sA'][::down2])
    pl.title('Primary Bridge Switching')
    pl.xlabel('t in (sec)')
    pl.ylabel("$s_{p}(t)$ (p.u)")
    pl.grid('on')

    pl.subplot(gs[1, 1])
    pl.plot(t[::down2], timeSw['sC'][::down2])
    pl.title('Secondary Bridge Switching')
    pl.xlabel('t in (sec)')
    pl.ylabel("$s_{s}(t)$ (p.u)")
    pl.grid('on')

    pl.subplot(gs[2, 0])
    pl.stem(f[::down][0:50], freqSw['Sa'][::down][0:50])
    pl.xlim(0, 50)
    pl.title('Primary Switching Spectrum')
    pl.xlabel("$f/f_{1}$ (Hz/Hz)")
    pl.ylabel("$S_{p}(f)$ (p.u)")
    pl.yscale('log')
    pl.ylim(0.0001, )
    pl.grid('on')

    pl.subplot(gs[2, 1])
    pl.stem(f[::down][0:50], freqSw['Xas'][::down][0:50])
    pl.xlim(0, 50)
    pl.title('Sampled Reference Spectrum')
    pl.xlabel("$f/f_{1}$ (Hz/Hz)")
    pl.ylabel("$X_{p}(f)$ (p.u)")
    pl.yscale('log')
    pl.ylim(0.0001, )
    pl.grid('on')

    ###################################################################################################################
    # Current
    ###################################################################################################################
    plt.figure()
    txt = "Currents: " + "$V_{dc}$=" + str(Vdc) + "V, " + "$M_{i}$=" + str(Mi) + "$ ,Q$=" + str(Q) + ", n=" + str(n_tr)
    plt.suptitle(txt, size=18)
    plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)

    i_pri = timeAc['i_ac_pri'] if 'i_ac_pri' in timeAc else timeAc['i_a']
    i_sec = timeAc['i_ac_sec'] if 'i_ac_sec' in timeAc else None
    plt.subplot(2, 3, 1)
    plt.plot(t[::down2], i_pri[::down2])
    if i_sec is not None:
        plt.plot(t[::down2], i_sec[::down2])
        plt.legend(["$i_{ac,pri}$", "$i_{ac,sec}$"], loc='upper right')
    plt.ylabel("$i_{ac}(t)$ (A)")
    plt.title('Time-domain AC Currents')
    plt.xlabel('time in (sec)')
    plt.grid('on')

    ax = plt.subplot(2, 3, 2)
    plt.stem(f[::down][0:50], freqAc['I_a'][::down][0:50])
    plt.ylabel("$I_{L}(f)$ (A)")
    plt.xlim(0, 50)
    plt.title('Frequency-domain Current')
    plt.xlabel("$f/f_{1}$ (Hz/Hz)")
    plt.yscale('log')
    plt.ylim(0.1 / OoM(max(freqAc['I_a'])), )
    txt = "THD=" + str(round(I_ph_thd * 100, 2)) + "%"
    plt.text(0.75, 0.90, txt, bbox=dict(facecolor='tab:blue', alpha=0.5), transform=ax.transAxes)
    plt.grid('on')

    if setup['Top']['sourceType'] == 'DAB':
        x_phi = M_i * 180 / np.pi
        x_lab = "Phase Shift (deg)"
    else:
        x_phi = M_i
        x_lab = "$M_{i}$ in (p.u)"

    plt.subplot(2, 3, 3)
    plt.plot(x_phi, distAc['num']['I_a_thd'])
    if setup['Exp']['plot'] == 2:
        plt.plot(x_phi, distAc['ana']['I_a_thd'], 'tab:blue', linestyle="", marker="o")
        plt.legend(['Numerical', 'Analytical'])
    plt.ylim(0, )
    plt.title('Distortion Current')
    plt.ylabel("$I_{L,rms}^{THD}$ (A)")
    plt.xlabel(x_lab)
    plt.grid('on')

    plt.subplot(2, 3, 4)
    plt.plot(t[::down2], timeDc['i_dc'][::down2], t[::down2],
             np.mean(timeDc['i_dc']) * np.ones(np.size(timeDc['i_dc'][::down2])), '--')
    plt.ylabel("$i_{dc}(t)$ (A)")
    plt.title('Time-domain DC Current')
    plt.xlabel('time in (sec)')
    plt.legend(["$i_{dc}$", "$I_{dc,avg}$"], loc='upper right')
    plt.grid('on')

    ax = plt.subplot(2, 3, 5)
    plt.stem(f[::down][0:50], freqDc['I_dc'][::down][0:50])
    plt.ylabel("$I_{dc}(f)$ (A)")
    plt.xlim(0, 50)
    plt.title('Frequency-domain DC Current')
    plt.xlabel("$f/f_{1}$ (Hz/Hz)")
    plt.yscale('log')
    plt.ylim(0.1 / OoM(max(freqDc['I_dc'])), )
    txt = "THD=" + str(round(I_dc_thd * 100, 2)) + "%"
    plt.text(0.75, 0.90, txt, bbox=dict(facecolor='tab:blue', alpha=0.5), transform=ax.transAxes)
    plt.grid('on')

    plt.subplot(2, 3, 6)
    plt.plot(x_phi, distDc['num']['I_dc_thd'])
    if setup['Exp']['plot'] == 2:
        plt.plot(x_phi, distDc['ana']['I_dc_thd'], 'tab:blue', linestyle="", marker="o")
        plt.legend(['Numerical', 'Analytical'])
    plt.ylim(0, )
    plt.title('Distortion DC Current')
    plt.ylabel("$I_{dc,rms}^{THD}$ (A)")
    plt.xlabel(x_lab)
    plt.grid('on')

    ###################################################################################################################
    # Voltage
    ###################################################################################################################
    plt.figure()
    txt = "Voltages: " + "$V_{dc}$=" + str(Vdc) + "V, " + "$M_{i}$=" + str(Mi) + ", n=" + str(n_tr)
    plt.suptitle(txt, size=18)
    plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)

    v_p = timeAc['v_a0']
    v_s = timeAc['v_b0']
    v_s_ref = timeAc['v_b0_ref'] if 'v_b0_ref' in timeAc else n_tr * v_s

    plt.subplot(2, 3, 1)
    plt.plot(t[::down2], timeAc['v_a'][::down2], t[::down2], v_p[::down2], t[::down2], v_s_ref[::down2])
    plt.ylabel("$v(t)$ (V)")
    plt.title('Time-domain DAB Voltages (Reflected)')
    plt.xlabel('time in (sec)')
    plt.legend(["$v_{p}-n v_{s}$", "$v_{p}$", "$n v_{s}$"], loc='upper right')
    plt.grid('on')

    ax = plt.subplot(2, 3, 2)
    plt.stem(f[::down][0:50], freqAc['V_a'][::down][0:50])
    plt.ylabel("$V(f)$ (V)")
    plt.xlim(0, 50)
    plt.title('Frequency-domain Voltage')
    plt.xlabel("$f/f_{1}$ (Hz/Hz)")
    plt.yscale('log')
    plt.ylim(0.1 / OoM(max(freqAc['V_a'])), )
    txt = "THD=" + str(round(V_ph_thd * 100, 2)) + "%"
    plt.text(0.75, 0.90, txt, bbox=dict(facecolor='tab:blue', alpha=0.5), transform=ax.transAxes)
    plt.grid('on')

    plt.subplot(2, 3, 3)
    plt.plot(x_phi, distAc['num']['V_a_thd'])
    if setup['Exp']['plot'] == 2:
        plt.plot(x_phi, distAc['ana']['V_a_thd'], 'tab:blue', linestyle="", marker="o")
        plt.legend(['Numerical', 'Analytical'])
    plt.title('Distortion Voltage')
    plt.ylabel("$V_{rms}^{THD}$ (V)")
    plt.xlabel(x_lab)
    plt.grid('on')

    plt.subplot(2, 3, 4)
    plt.plot(t[::down2], timeDc['v_in'][::down2], t[::down2], timeDc['v_dc'][::down2], t[::down2],
             np.mean(timeDc['v_dc']) * np.ones(np.size(timeDc['v_dc'][::down2])), '--')
    plt.ylabel("$v_{dc}(t)$ (V)")
    plt.title('Time-domain DC Voltage')
    plt.xlabel('time in (sec)')
    plt.legend(["$v_{in}$", "$v_{dc}$", "$V_{dc,avg}$"], loc='upper right')
    plt.grid('on')

    ax = plt.subplot(2, 3, 5)
    plt.stem(f[::down][0:50], freqDc['V_dc'][::down][0:50])
    plt.ylabel("$V_{dc}(f)$ (A)")
    plt.xlim(0, 50)
    plt.title('Frequency-domain DC Voltage')
    plt.xlabel("$f/f_{1}$ (Hz/Hz)")
    plt.yscale('log')
    plt.ylim(0.1 / OoM(max(freqDc['V_dc'])), )
    txt = "THD=" + str(round(V_dc_thd * 100, 2)) + "%"
    plt.text(0.75, 0.90, txt, bbox=dict(facecolor='tab:blue', alpha=0.5), transform=ax.transAxes)
    plt.grid('on')

    plt.subplot(2, 3, 6)
    plt.plot(x_phi, distDc['num']['V_dc_thd'])
    if setup['Exp']['plot'] == 2:
        plt.plot(x_phi, distDc['ana']['V_dc_thd'], 'tab:blue', linestyle="", marker="o")
        plt.legend(['Numerical', 'Analytical'])
    plt.title('Distortion DC Voltage')
    plt.ylabel("$V_{dc,rms}^{THD}$ (V)")
    plt.xlabel(x_lab)
    plt.grid('on')

    ###################################################################################################################
    # Power vs Phase Shift (Mi)
    ###################################################################################################################
    plt.figure()
    if setup['Top']['sourceType'] == 'DAB':
        phi_deg = M_i * 180 / np.pi
    else:
        phi_deg = M_i * 90
    S_app = distAc['num']['V_a_eff'] * distAc['num']['I_a_eff']
    plt.plot(phi_deg, S_app)
    plt.title('Apparent Power vs Phase Shift')
    plt.xlabel("Phase Shift (deg)")
    plt.ylabel("Apparent Power (VA)")
    plt.grid('on')
    plt.show()

    ###################################################################################################################
    # MSG OUT
    ###################################################################################################################
    print("END: Plotting DAB Sweep Results")
