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
    txt = "DAB Switching Strategies: " + "$M_{i}$=" + str(Mi) + "$ ,Q$=" + str(Q) + ", Ntr=" + str(n_tr)
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
    pl.xlabel("$f/f_{1}$")
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
    pl.xlabel("$f/f_{1}$")
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
    # 4) Distortion AC Side
    ###################################################################################################################
    plt.figure(figsize=(12, 8))
    plt.suptitle("DAB Distortion Overview: AC Side", size=18)
    plt.subplots_adjust(hspace=0.4, wspace=0.3)

    if setup['Top']['sourceType'] == 'DAB':
        x_phi = M_i * 180 / np.pi
        x_lab = "Phase Shift (deg)"
    else:
        x_phi = M_i
        x_lab = "$M_{i}$ in (p.u)"

    # Voltage THD
    plt.subplot(2, 2, 1)
    plt.plot(x_phi, distAc['num']['V_a_thd'] * 100)
    plt.title('Primary Voltage THD')
    plt.ylabel('THD (%)')
    plt.xlabel(x_lab)
    plt.grid('on')

    # Current THD
    plt.subplot(2, 2, 2)
    plt.plot(x_phi, distAc['num']['I_a_thd'] * 100)
    plt.title('Primary Current THD')
    plt.ylabel('THD (%)')
    plt.xlabel(x_lab)
    plt.grid('on')

    # Apparent Power
    plt.subplot(2, 2, 3)
    S_app = distAc['num']['V_a_eff'] * distAc['num']['I_a_eff']
    plt.plot(x_phi, S_app)
    plt.title('Apparent Power vs Phase Shift')
    plt.xlabel(x_lab)
    plt.ylabel("Apparent Power (VA)")
    plt.grid('on')

    # Power Factor
    plt.subplot(2, 2, 4)
    # Placeholder for PF if not explicitly in distAc
    plt.plot(x_phi, np.ones(len(x_phi)))
    plt.title('Power Factor (Placeholder)')
    plt.xlabel(x_lab)
    plt.ylabel("PF (-)")
    plt.grid('on')

    ###################################################################################################################
    # 5) Distortion DC Side
    ###################################################################################################################
    plt.figure(figsize=(12, 8))
    plt.suptitle("DAB Distortion Overview: DC Side", size=18)
    plt.subplots_adjust(hspace=0.4, wspace=0.3)

    # Voltage THD
    plt.subplot(2, 2, 1)
    plt.plot(x_phi, distDc['num']['V_dc_thd'] * 100)
    plt.title('DC Voltage THD')
    plt.ylabel('THD (%)')
    plt.xlabel(x_lab)
    plt.grid('on')

    # Current THD
    plt.subplot(2, 2, 2)
    plt.plot(x_phi, distDc['num']['I_dc_thd'] * 100)
    plt.title('DC Current THD')
    plt.ylabel('THD (%)')
    plt.xlabel(x_lab)
    plt.grid('on')

    # DC Power
    plt.subplot(2, 2, 3)
    P_dc = distDc['num']['V_dc_eff'] * distDc['num']['I_dc_eff']
    plt.plot(x_phi, P_dc)
    plt.title('DC Power vs Phase Shift')
    plt.xlabel(x_lab)
    plt.ylabel("Power (W)")
    plt.grid('on')
    plt.show()

    ###################################################################################################################
    # MSG OUT
    ###################################################################################################################
    print("END: Plotting DAB Sweep Results")
