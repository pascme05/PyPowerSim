#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         plotStat_DAB
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
def plotStat_DAB(time, freq, setup):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("START: Plotting Stationary DAB Waveforms")

    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    timeSw = time['Sw']
    timeAc = time['Ac']
    timeDc = time['Dc']
    timeElec = time.get('Elec', {})
    timeLoss = time.get('Loss', {})
    timeTher = time.get('Ther', {})
    freqSw = freq['Sw']
    freqAc = freq['Ac']
    freqDc = freq['Dc']

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
    # Load angle
    ###################################################################################################################
    Y = fft(timeAc['v_a'])
    phiV = np.angle(Y)[1]
    Y = fft(timeAc['i_a'])
    phiI = np.angle(Y)[1]
    phi = phiV - phiI + 2 * np.pi
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
    txt = "DAB Switching Functions: " + "$M_{i}$=" + str(Mi) + "$ ,Q$=" + str(Q) + ", n=" + str(n_tr)
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
    pl.plot(t[::down2], timeSw['sA'][::down2], label='S_A')
    pl.plot(t[::down2], timeSw['sB'][::down2], label='S_B')
    pl.title('Primary Bridge Switching')
    pl.xlabel('t in (sec)')
    pl.ylabel("$s(t)$ (p.u)")
    pl.legend(loc='upper right')
    pl.grid('on')

    pl.subplot(gs[1, 1])
    pl.plot(t[::down2], timeSw['sC'][::down2], label='S_C')
    pl.plot(t[::down2], timeSw['sD'][::down2], label='S_D')
    pl.title('Secondary Bridge Switching')
    pl.xlabel('t in (sec)')
    pl.ylabel("$s(t)$ (p.u)")
    pl.legend(loc='upper right')
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
    # Currents
    ###################################################################################################################
    plt.figure()
    txt = "DAB Currents: " + "$V_{dc}$=" + str(Vdc) + "V, " + "$M_{i}$=" + str(Mi) + ", n=" + str(n_tr)
    plt.suptitle(txt, size=18)
    plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)

    plt.subplot(2, 2, 1)
    i_pri = timeAc['i_ac_pri'] if 'i_ac_pri' in timeAc else timeAc['i_a']
    i_sec = timeAc['i_ac_sec'] if 'i_ac_sec' in timeAc else None
    plt.plot(t[::down2], i_pri[::down2])
    if i_sec is not None:
        plt.plot(t[::down2], i_sec[::down2])
        plt.legend(["$i_{ac,pri}$", "$i_{ac,sec}$"], loc='upper right')
    plt.ylabel("$i_{ac}(t)$ (A)")
    plt.title('Time-domain AC Currents')
    plt.xlabel('time in (sec)')
    plt.grid('on')

    ax = plt.subplot(2, 2, 2)
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

    plt.subplot(2, 2, 3)
    plt.plot(t[::down2], timeDc['i_dc'][::down2], t[::down2],
             np.mean(timeDc['i_dc']) * np.ones(np.size(timeDc['i_dc'][::down2])), '--')
    plt.ylabel("$i_{dc}(t)$ (A)")
    plt.title('Time-domain DC Current')
    plt.xlabel('time in (sec)')
    plt.legend(["$i_{dc}$", "$I_{dc,avg}$"], loc='upper right')
    plt.grid('on')

    ax = plt.subplot(2, 2, 4)
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

    ###################################################################################################################
    # Voltages
    ###################################################################################################################
    plt.figure()
    txt = "DAB Voltages: " + "$V_{dc}$=" + str(Vdc) + "V, " + "$M_{i}$=" + str(Mi) + ", n=" + str(n_tr)
    plt.suptitle(txt, size=18)
    plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)

    plt.subplot(2, 2, 1)
    v_p = timeAc['v_a0']
    v_s = timeAc['v_b0']
    v_s_ref = timeAc['v_b0_ref'] if 'v_b0_ref' in timeAc else n_tr * v_s
    plt.plot(t[::down2], timeAc['v_a'][::down2], t[::down2], v_p[::down2], t[::down2], v_s_ref[::down2])
    plt.ylabel("$v(t)$ (V)")
    plt.title('Time-domain DAB Voltages (Reflected)')
    plt.xlabel('time in (sec)')
    plt.legend(["$v_{p}-n v_{s}$", "$v_{p}$", "$n v_{s}$"], loc='upper right')
    plt.grid('on')

    ax = plt.subplot(2, 2, 2)
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

    plt.subplot(2, 2, 3)
    plt.plot(t[::down2], timeDc['v_in'][::down2], t[::down2], timeDc['v_dc'][::down2], t[::down2],
             np.mean(timeDc['v_dc']) * np.ones(np.size(timeDc['v_dc'][::down2])), '--')
    plt.ylabel("$v_{dc}(t)$ (V)")
    plt.title('Time-domain DC Voltage')
    plt.xlabel('time in (sec)')
    plt.legend(["$v_{in}$", "$v_{dc}$", "$V_{dc,avg}$"], loc='upper right')
    plt.grid('on')

    ax = plt.subplot(2, 2, 4)
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

    ###################################################################################################################
    # Losses and Thermal (Switches)
    ###################################################################################################################
    if 'sw' in timeElec and len(timeElec['sw']) > 0:
        def _key(x):
            try:
                return int(x[1:])
            except:
                return 0

        def _x_from(y):
            y_len = len(y)
            if y_len <= 1:
                return np.array([0])
            return np.linspace(t[0], t[-1], y_len)

        sw_ids = sorted(list(timeElec['sw'].keys()), key=_key)
        s_primary = 'S1' if 'S1' in sw_ids else sw_ids[0]
        s_secondary = 'S5' if 'S5' in sw_ids else sw_ids[-1]

        plt.figure()
        txt = "DAB Switch Losses (Primary/Secondary): " + s_primary + " & " + s_secondary
        plt.suptitle(txt, size=18)
        plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)

        # Currents
        x_i_tp = _x_from(timeElec['sw'][s_primary]['i_T'])
        x_i_dp = _x_from(timeElec['sw'][s_primary]['i_D'])
        x_i_ts = _x_from(timeElec['sw'][s_secondary]['i_T'])
        x_i_ds = _x_from(timeElec['sw'][s_secondary]['i_D'])

        plt.subplot(4, 2, 1)
        plt.plot(x_i_tp[::down2], timeElec['sw'][s_primary]['i_T'][::down2], x_i_dp[::down2], timeElec['sw'][s_primary]['i_D'][::down2])
        plt.ylabel("$i(t)$ (A)")
        plt.title('Currents ' + s_primary)
        plt.xticks([], [])
        plt.legend(["$i_{T}$", "$i_{D}$"], loc='upper right')
        plt.grid('on')

        plt.subplot(4, 2, 2)
        plt.plot(x_i_ts[::down2], timeElec['sw'][s_secondary]['i_T'][::down2], x_i_ds[::down2], timeElec['sw'][s_secondary]['i_D'][::down2])
        plt.ylabel("$i(t)$ (A)")
        plt.title('Currents ' + s_secondary)
        plt.xticks([], [])
        plt.legend(["$i_{T}$", "$i_{D}$"], loc='upper right')
        plt.grid('on')

        # Voltages
        x_v_tp = _x_from(timeElec['sw'][s_primary]['v_T'])
        x_v_dp = _x_from(timeElec['sw'][s_primary]['v_D'])
        x_v_ts = _x_from(timeElec['sw'][s_secondary]['v_T'])
        x_v_ds = _x_from(timeElec['sw'][s_secondary]['v_D'])

        plt.subplot(4, 2, 3)
        plt.plot(x_v_tp[::down2], timeElec['sw'][s_primary]['v_T'][::down2], x_v_dp[::down2], timeElec['sw'][s_primary]['v_D'][::down2])
        plt.ylabel("$v(t)$ (V)")
        plt.title('Voltages ' + s_primary)
        plt.xticks([], [])
        plt.legend(["$v_{T}$", "$v_{D}$"], loc='upper right')
        plt.grid('on')

        plt.subplot(4, 2, 4)
        plt.plot(x_v_ts[::down2], timeElec['sw'][s_secondary]['v_T'][::down2], x_v_ds[::down2], timeElec['sw'][s_secondary]['v_D'][::down2])
        plt.ylabel("$v(t)$ (V)")
        plt.title('Voltages ' + s_secondary)
        plt.xticks([], [])
        plt.legend(["$v_{T}$", "$v_{D}$"], loc='upper right')
        plt.grid('on')

        # Conduction losses
        if 'sw' in timeLoss:
            x_p_tp = _x_from(timeLoss['sw'][s_primary]['p_T_c'])
            x_p_dp = _x_from(timeLoss['sw'][s_primary]['p_D_c'])
            x_p_ts = _x_from(timeLoss['sw'][s_secondary]['p_T_c'])
            x_p_ds = _x_from(timeLoss['sw'][s_secondary]['p_D_c'])

            plt.subplot(4, 2, 5)
            plt.plot(x_p_tp[::down2], timeLoss['sw'][s_primary]['p_T_c'][::down2], x_p_dp[::down2], timeLoss['sw'][s_primary]['p_D_c'][::down2])
            plt.ylabel("$p(t)$ (W)")
            plt.title('Conduction Losses ' + s_primary)
            plt.xticks([], [])
            plt.legend(["$p_{T}$", "$p_{D}$"], loc='upper right')
            plt.grid('on')

            plt.subplot(4, 2, 6)
            plt.plot(x_p_ts[::down2], timeLoss['sw'][s_secondary]['p_T_c'][::down2], x_p_ds[::down2], timeLoss['sw'][s_secondary]['p_D_c'][::down2])
            plt.ylabel("$p(t)$ (W)")
            plt.title('Conduction Losses ' + s_secondary)
            plt.xticks([], [])
            plt.legend(["$p_{T}$", "$p_{D}$"], loc='upper right')
            plt.grid('on')

            # Switching losses
            x_ps_tp = _x_from(timeLoss['sw'][s_primary]['p_T_s'])
            x_ps_dp = _x_from(timeLoss['sw'][s_primary]['p_D_s'])
            x_ps_ts = _x_from(timeLoss['sw'][s_secondary]['p_T_s'])
            x_ps_ds = _x_from(timeLoss['sw'][s_secondary]['p_D_s'])

            plt.subplot(4, 2, 7)
            plt.plot(x_ps_tp[::down2], timeLoss['sw'][s_primary]['p_T_s'][::down2], x_ps_dp[::down2], timeLoss['sw'][s_primary]['p_D_s'][::down2])
            plt.ylabel("$p(t)$ (W)")
            plt.title('Switching Losses ' + s_primary)
            plt.xlabel('time in (sec)')
            plt.legend(["$p_{T}$", "$p_{D}$"], loc='upper right')
            plt.grid('on')

            plt.subplot(4, 2, 8)
            plt.plot(x_ps_ts[::down2], timeLoss['sw'][s_secondary]['p_T_s'][::down2], x_ps_ds[::down2], timeLoss['sw'][s_secondary]['p_D_s'][::down2])
            plt.ylabel("$p(t)$ (W)")
            plt.title('Switching Losses ' + s_secondary)
            plt.xlabel('time in (sec)')
            plt.legend(["$p_{T}$", "$p_{D}$"], loc='upper right')
            plt.grid('on')

        # Thermal
        if 'sw' in timeTher:
            def _pick(ids, want):
                return want if want in ids else ids[0]

            s_p_hs = _pick(sw_ids, 'S1')
            s_p_ls = _pick(sw_ids, 'S2') if len(sw_ids) > 1 else s_p_hs
            s_s_hs = _pick(sw_ids, 'S5') if len(sw_ids) > 4 else sw_ids[-2]
            s_s_ls = _pick(sw_ids, 'S6') if len(sw_ids) > 5 else sw_ids[-1]

            def _therm_axes(sid):
                x_t = _x_from(timeTher['sw']['T' + sid[1:]])
                x_d = _x_from(timeTher['sw']['D' + sid[1:]])
                x_c = _x_from(timeTher['sw']['C' + sid[1:]])
                return x_t, x_d, x_c

            plt.figure()
            plt.suptitle("DAB Thermal (Primary/Secondary, HS/LS)", size=18)
            plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)

            # Primary HS
            x_t, x_d, x_c = _therm_axes(s_p_hs)
            plt.subplot(2, 2, 1)
            plt.plot(x_t[::down2], timeTher['sw']['T' + s_p_hs[1:]][::down2],
                     x_d[::down2], timeTher['sw']['D' + s_p_hs[1:]][::down2],
                     x_c[::down2], timeTher['sw']['C' + s_p_hs[1:]][::down2])
            plt.ylabel("$\\Theta(t)$ (C)")
            plt.title('Primary HS ' + s_p_hs)
            plt.grid('on')

            # Primary LS
            x_t, x_d, x_c = _therm_axes(s_p_ls)
            plt.subplot(2, 2, 2)
            plt.plot(x_t[::down2], timeTher['sw']['T' + s_p_ls[1:]][::down2],
                     x_d[::down2], timeTher['sw']['D' + s_p_ls[1:]][::down2],
                     x_c[::down2], timeTher['sw']['C' + s_p_ls[1:]][::down2])
            plt.ylabel("$\\Theta(t)$ (C)")
            plt.title('Primary LS ' + s_p_ls)
            plt.grid('on')

            # Secondary HS
            x_t, x_d, x_c = _therm_axes(s_s_hs)
            plt.subplot(2, 2, 3)
            plt.plot(x_t[::down2], timeTher['sw']['T' + s_s_hs[1:]][::down2],
                     x_d[::down2], timeTher['sw']['D' + s_s_hs[1:]][::down2],
                     x_c[::down2], timeTher['sw']['C' + s_s_hs[1:]][::down2])
            plt.ylabel("$\\Theta(t)$ (C)")
            plt.title('Secondary HS ' + s_s_hs)
            plt.grid('on')

            # Secondary LS
            x_t, x_d, x_c = _therm_axes(s_s_ls)
            plt.subplot(2, 2, 4)
            plt.plot(x_t[::down2], timeTher['sw']['T' + s_s_ls[1:]][::down2],
                     x_d[::down2], timeTher['sw']['D' + s_s_ls[1:]][::down2],
                     x_c[::down2], timeTher['sw']['C' + s_s_ls[1:]][::down2])
            plt.ylabel("$\\Theta(t)$ (C)")
            plt.title('Secondary LS ' + s_s_ls)
            plt.grid('on')

    plt.show()

    ###################################################################################################################
    # MSG OUT
    ###################################################################################################################
    print("END: Plotting Stationary DAB Waveforms")
