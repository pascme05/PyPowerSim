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

    ###################################################################################################################
    # Switching Functions
    ###################################################################################################################
    gs = gridspec.GridSpec(3, 2)
    pl.figure()
    txt = "DAB Switching (Transient): " + "$M_{i}$=" + str(Mi) + "$ ,Q$=" + str(Q) + ", n=" + str(n_tr)
    pl.suptitle(txt, size=18)
    pl.subplots_adjust(hspace=0.35, wspace=0.20, left=0.075, right=0.925, top=0.90, bottom=0.075)

    ax1 = pl.subplot(gs[0, :])
    ax2 = ax1.twinx()
    ax1.plot(timeSw['t'][::down2], timeSw['c'][::down2], 'black')
    ax2.plot(timeSw['t'][::down2], timeSw['sA'][::down2], color='tab:blue')
    ax2.plot(timeSw['t'][::down2], timeSw['sC'][::down2], color='tab:orange')
    pl.title('Primary and Secondary Switching States')
    pl.xlabel('t in (sec)')
    ax1.set_ylabel('c(t)', color='black')
    ax2.set_ylabel('s(t) (p.u)', color='tab:blue')
    pl.legend(["$c$", "$s_{p}$", "$s_{s}$"], loc='upper right')
    pl.grid('on')

    pl.subplot(gs[1, 0])
    pl.plot(timeSw['t'][::down2], timeSw['sA'][::down2])
    pl.title('Primary Bridge Switching')
    pl.xlabel('t in (sec)')
    pl.ylabel("$s_{p}(t)$ (p.u)")
    pl.grid('on')

    pl.subplot(gs[1, 1])
    pl.plot(timeSw['t'][::down2], timeSw['sC'][::down2])
    pl.title('Secondary Bridge Switching')
    pl.xlabel('t in (sec)')
    pl.ylabel("$s_{s}(t)$ (p.u)")
    pl.grid('on')

    pl.subplot(gs[2, 0])
    n1 = min(len(f[::5]), len(freqSw['Sa'][::5]), 50)
    pl.stem(f[::5][0:n1], freqSw['Sa'][::5][0:n1])
    pl.xlim(0, 50)
    pl.title('Primary Switching Spectrum')
    pl.xlabel("$f/f_{1}$ (Hz/Hz)")
    pl.ylabel("$S_{p}(f)$ (p.u)")
    pl.yscale('log')
    pl.ylim(0.0001, )
    pl.grid('on')

    pl.subplot(gs[2, 1])
    n2 = min(len(f[::5]), len(freqSw['Xas'][::5]), 50)
    pl.stem(f[::5][0:n2], freqSw['Xas'][::5][0:n2])
    pl.xlim(0, 50)
    pl.title('Sampled Reference Spectrum')
    pl.xlabel("$f/f_{1}$ (Hz/Hz)")
    pl.ylabel("$X_{p}(f)$ (p.u)")
    pl.yscale('log')
    pl.ylim(0.0001, )
    pl.grid('on')

    ###################################################################################################################
    # Currents and Voltages
    ###################################################################################################################
    plt.figure()
    txt = "DAB Transient Currents/Voltages: " + "$V_{dc}$=" + str(Vdc) + "V, n=" + str(n_tr)
    plt.suptitle(txt, size=18)
    plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)

    i_pri = timeAc['i_ac_pri'] if 'i_ac_pri' in timeAc else timeAc['i_a']
    i_sec = timeAc['i_ac_sec'] if 'i_ac_sec' in timeAc else None
    x_i = _x_from(i_pri)
    plt.subplot(2, 2, 1)
    plt.plot(x_i[::down2], i_pri[::down2])
    if i_sec is not None:
        x_i2 = _x_from(i_sec)
        plt.plot(x_i2[::down2], i_sec[::down2])
        plt.legend(["$i_{ac,pri}$", "$i_{ac,sec}$"], loc='upper right')
    plt.ylabel("$i_{ac}(t)$ (A)")
    plt.title('AC Currents')
    plt.xlabel('time in (sec)')
    plt.grid('on')

    plt.subplot(2, 2, 2)
    v_p = timeAc['v_a0']
    v_s = timeAc['v_b0']
    v_s_ref = timeAc['v_b0_ref'] if 'v_b0_ref' in timeAc else n_tr * v_s
    x_v = _x_from(timeAc['v_a'])
    x_vp = _x_from(v_p)
    x_vs = _x_from(v_s_ref)
    plt.plot(x_v[::down2], timeAc['v_a'][::down2], x_vp[::down2], v_p[::down2], x_vs[::down2], v_s_ref[::down2])
    plt.ylabel("$v(t)$ (V)")
    plt.title('DAB Voltages (Reflected)')
    plt.xlabel('time in (sec)')
    plt.legend(["$v_{p}-n v_{s}$", "$v_{p}$", "$n v_{s}$"], loc='upper right')
    plt.grid('on')

    x_idc = _x_from(timeDc['i_dc'])
    plt.subplot(2, 2, 3)
    plt.plot(x_idc[::down2], timeDc['i_dc'][::down2])
    plt.ylabel("$i_{dc}(t)$ (A)")
    plt.title('DC Current')
    plt.xlabel('time in (sec)')
    plt.grid('on')

    x_vdc = _x_from(timeDc['v_dc'])
    plt.subplot(2, 2, 4)
    plt.plot(x_vdc[::down2], timeDc['v_dc'][::down2])
    plt.ylabel("$v_{dc}(t)$ (V)")
    plt.title('DC Voltage')
    plt.xlabel('time in (sec)')
    plt.grid('on')

    ###################################################################################################################
    # Losses (Primary/Secondary, HS/LS)
    ###################################################################################################################
    if 'sw' in timeLoss and len(timeLoss['sw']) > 0:
        def _key(x):
            try:
                return int(x[1:])
            except:
                return 0

        sw_ids = sorted(list(timeLoss['sw'].keys()), key=_key)

        def _pick(ids, want):
            return want if want in ids else ids[0]

        s_p_hs = _pick(sw_ids, 'S1')
        s_p_ls = _pick(sw_ids, 'S2') if len(sw_ids) > 1 else s_p_hs
        s_s_hs = _pick(sw_ids, 'S5') if len(sw_ids) > 4 else sw_ids[-2]
        s_s_ls = _pick(sw_ids, 'S6') if len(sw_ids) > 5 else sw_ids[-1]

        def _loss_axes(sid):
            return np.linspace(t[0], t[-1], len(timeLoss['sw'][sid]['p_T']))

        plt.figure()
        plt.suptitle("DAB Losses (Transient)", size=18)
        plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)

        # Primary HS
        x_t = _loss_axes(s_p_hs)
        plt.subplot(2, 2, 1)
        plt.plot(x_t[::down2], timeLoss['sw'][s_p_hs]['p_T'][::down2],
                 x_t[::down2], timeLoss['sw'][s_p_hs]['p_D'][::down2])
        plt.ylabel("$p(t)$ (W)")
        plt.title('Primary HS ' + s_p_hs)
        plt.legend(["$p_{T}$", "$p_{D}$"], loc='upper right')
        plt.grid('on')

        # Primary LS
        x_t = _loss_axes(s_p_ls)
        plt.subplot(2, 2, 2)
        plt.plot(x_t[::down2], timeLoss['sw'][s_p_ls]['p_T'][::down2],
                 x_t[::down2], timeLoss['sw'][s_p_ls]['p_D'][::down2])
        plt.ylabel("$p(t)$ (W)")
        plt.title('Primary LS ' + s_p_ls)
        plt.legend(["$p_{T}$", "$p_{D}$"], loc='upper right')
        plt.grid('on')

        # Secondary HS
        x_t = _loss_axes(s_s_hs)
        plt.subplot(2, 2, 3)
        plt.plot(x_t[::down2], timeLoss['sw'][s_s_hs]['p_T'][::down2],
                 x_t[::down2], timeLoss['sw'][s_s_hs]['p_D'][::down2])
        plt.ylabel("$p(t)$ (W)")
        plt.title('Secondary HS ' + s_s_hs)
        plt.legend(["$p_{T}$", "$p_{D}$"], loc='upper right')
        plt.grid('on')

        # Secondary LS
        x_t = _loss_axes(s_s_ls)
        plt.subplot(2, 2, 4)
        plt.plot(x_t[::down2], timeLoss['sw'][s_s_ls]['p_T'][::down2],
                 x_t[::down2], timeLoss['sw'][s_s_ls]['p_D'][::down2])
        plt.ylabel("$p(t)$ (W)")
        plt.title('Secondary LS ' + s_s_ls)
        plt.legend(["$p_{T}$", "$p_{D}$"], loc='upper right')
        plt.grid('on')

    ###################################################################################################################
    # Thermal (Primary/Secondary, HS/LS)
    ###################################################################################################################
    if 'sw' in timeTher and len(timeTher['sw']) > 0:
        def _key(x):
            try:
                return int(x[1:])
            except:
                return 0

        sw_ids = sorted(list(timeTher['sw'].keys()), key=_key)

        def _pick(ids, want):
            return want if want in ids else ids[0]

        s_p_hs = _pick(sw_ids, 'T1')
        s_p_ls = _pick(sw_ids, 'T2') if len(sw_ids) > 1 else s_p_hs
        s_s_hs = _pick(sw_ids, 'T5') if len(sw_ids) > 4 else sw_ids[-2]
        s_s_ls = _pick(sw_ids, 'T6') if len(sw_ids) > 5 else sw_ids[-1]

        def _therm_axes(sid):
            return np.linspace(t[0], t[-1], len(timeTher['sw'][sid]))

        plt.figure()
        plt.suptitle("DAB Thermal (Transient)", size=18)
        plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)

        # Primary HS
        x_t = _therm_axes(s_p_hs)
        plt.subplot(2, 2, 1)
        plt.plot(x_t[::down2], timeTher['sw'][s_p_hs][::down2])
        plt.ylabel("$\\Theta(t)$ (C)")
        plt.title('Primary HS ' + s_p_hs)
        plt.grid('on')

        # Primary LS
        x_t = _therm_axes(s_p_ls)
        plt.subplot(2, 2, 2)
        plt.plot(x_t[::down2], timeTher['sw'][s_p_ls][::down2])
        plt.ylabel("$\\Theta(t)$ (C)")
        plt.title('Primary LS ' + s_p_ls)
        plt.grid('on')

        # Secondary HS
        x_t = _therm_axes(s_s_hs)
        plt.subplot(2, 2, 3)
        plt.plot(x_t[::down2], timeTher['sw'][s_s_hs][::down2])
        plt.ylabel("$\\Theta(t)$ (C)")
        plt.title('Secondary HS ' + s_s_hs)
        plt.grid('on')

        # Secondary LS
        x_t = _therm_axes(s_s_ls)
        plt.subplot(2, 2, 4)
        plt.plot(x_t[::down2], timeTher['sw'][s_s_ls][::down2])
        plt.ylabel("$\\Theta(t)$ (C)")
        plt.title('Secondary LS ' + s_s_ls)
        plt.grid('on')

    plt.show()

    ###################################################################################################################
    # MSG OUT
    ###################################################################################################################
    print("END: Plotting Transient DAB Waveforms")
