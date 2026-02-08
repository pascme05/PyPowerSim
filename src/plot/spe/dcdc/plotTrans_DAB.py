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

    timeSw = time['Sw']
    timeAc = time['Ac']
    timeDc = time['Dc']
    timeLoss = time.get('Loss', {})
    timeTher = time.get('Ther', {})
    freqSw = freq['Sw']
    freqAc = freq['Ac']
    freqDc = freq['Dc']

    fel = setup['Top']['fel']
    fsim = setup['Exp']['fsim']
    fs = setup['Par']['PWM']['fs']
    Q = int(fs / fel)
    Mi = setup['Dat']['stat']['Mi']
    Vdc = setup['Dat']['stat']['Vdc']
    n_tr = setup['Top'].get('n', 1)
    down = int(setup['Dat']['stat']['cyc']) - 2
    if down < 1:
        down = 1

    if 't' in time:
        t = time['t']
    else:
        t = timeSw['t'].values

    def _x_from(y):
        y_len = len(y)
        if y_len <= 1:
            return np.array([0])
        return np.linspace(t[0], t[-1], y_len)

    down2 = int(fsim / fs / 200)
    if down2 < 1:
        down2 = 1

    f = fsim * np.linspace(0, 0.5, int(len(timeSw['t']) / 2)) / fel

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
    txt = "DAB Switching Strategies (Transient): " + "$M_{i}$=" + str(Mi) + "$ ,Q$=" + str(Q) + ", Ntr=" + str(n_tr)
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
    n1 = min(len(f_plot[::5]), len(freqSw['Sa'][::5]), 50)
    pl.stem(f_plot[::5][0:n1], freqSw['Sa'][::5][0:n1])
    pl.xlim(0, 50)
    pl.title('Primary Switching Spectrum')
    pl.xlabel("$f/f_{1}$")
    pl.ylabel("$S_{p}(f)$ (p.u)")
    pl.yscale('log')
    pl.ylim(0.0001, )
    pl.grid('on')

    # Secondary
    pl.subplot(gs[2, 1])
    freq_sec = freqSw.get('Sb', freqSw.get('Sc'))
    if freq_sec is not None:
        n2 = min(len(f_plot[::5]), len(freq_sec[::5]), 50)
        pl.stem(f_plot[::5][0:n2], freq_sec[::5][0:n2])
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
    plt.suptitle("DAB Primary Side: Current and Voltage (Transient)", size=18)
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
    plt.ylim(0.0001, )
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
    plt.ylim(0.0001, )
    plt.grid('on')

    ###################################################################################################################
    # 3) Secondary Side Current and Voltage
    ###################################################################################################################
    plt.figure(figsize=(12, 8))
    plt.suptitle("DAB Secondary Side: Current and Voltage (Transient)", size=18)
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
    plt.ylim(0.0001, )
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
    plt.ylim(0.0001, )
    plt.grid('on')

    ###################################################################################################################
    # 4) Power Losses Overview (Transient)
    ###################################################################################################################
    if 'sw' in timeLoss and len(timeLoss['sw']) > 0:
        plt.figure(figsize=(12, 8))
        plt.suptitle("DAB Losses Overview (Transient)", size=18)
        plt.subplots_adjust(hspace=0.4, wspace=0.3)

        def _key(x):
            try: return int(x[1:])
            except: return 0
        sw_ids = sorted(list(timeLoss['sw'].keys()), key=_key)

        # Primary HS and LS
        plt.subplot(2, 2, 1)
        for sid in sw_ids[:2]:
            plt.plot(t_plot[::down2], timeLoss['sw'][sid]['p_L'][::down2], label=sid)
        plt.title('Primary Bridge Losses')
        plt.xlabel(t_label)
        plt.ylabel('Power (W)')
        plt.legend()
        plt.grid('on')

        # Secondary HS and LS
        plt.subplot(2, 2, 2)
        for sid in sw_ids[4:6]:
            plt.plot(t_plot[::down2], timeLoss['sw'][sid]['p_L'][::down2], label=sid)
        plt.title('Secondary Bridge Losses')
        plt.xlabel(t_label)
        plt.ylabel('Power (W)')
        plt.legend()
        plt.grid('on')

        # Total Losses
        plt.subplot(2, 1, 2)
        p_total = np.zeros(len(t_plot))
        for sid in timeLoss['sw']:
            p_total += timeLoss['sw'][sid]['p_L'].values
        plt.plot(t_plot[::down2], p_total[::down2])
        plt.title('Total Switch Losses')
        plt.xlabel(t_label)
        plt.ylabel('Power (W)')
        plt.grid('on')

    ###################################################################################################################
    # 5) Thermal Overview (Transient)
    ###################################################################################################################
    if 'sw' in timeTher and len(timeTher['sw']) > 0:
        plt.figure(figsize=(12, 8))
        plt.suptitle("DAB Thermal Overview (Transient)", size=18)
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

        # Capacitor
        plt.subplot(2, 2, 3)
        if 'cap' in timeTher:
            for cid in timeTher['cap']:
                plt.plot(t_plot[::down2], timeTher['cap'][cid][::down2], label=cid)
        plt.title('Capacitor Temperature')
        plt.xlabel(t_label)
        plt.ylabel('Temperature (°C)')
        plt.legend()
        plt.grid('on')

    plt.show()

    ###################################################################################################################
    # MSG OUT
    ###################################################################################################################
    print("END: Plotting Transient DAB Waveforms")
