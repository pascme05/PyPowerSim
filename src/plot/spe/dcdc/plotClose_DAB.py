#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         plotClose_DAB
# Date:         14.02.2026
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
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')


#######################################################################################################################
# Function
#######################################################################################################################
def plotClose_DAB(time, freq, setup):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("START: Plotting Closed Loop DAB Waveforms")

    # ==============================================================================
    # Init
    # ==============================================================================
    timeSw = time.get('Sw', {})
    timeAc = time.get('Ac', {})
    timeDc = time.get('Dc', {})

    # Time base
    if 't' in time:
        t_main = time['t']
    elif isinstance(timeSw, dict) and 't' in timeSw:
        t_main = timeSw['t'].values
    else:
        t_main = np.linspace(0, setup['Dat']['trans']['tmax'], len(next(iter(timeDc.values()))) if timeDc else 1)

    def _time_for(sig):
        n = len(sig)
        if n == len(t_main):
            return t_main
        return np.linspace(t_main[0], t_main[-1], n)

    # Downsample
    fs = setup['Par']['PWM']['fs']
    fsim = setup['Exp']['fsim']
    down2 = int(fsim / fs / 200)
    if down2 < 1:
        down2 = 1

    ###################################################################################################################
    # Reference tracking
    ###################################################################################################################
    i_ref = None
    if isinstance(timeDc, dict) and 'i_ref' in timeDc:
        i_ref = timeDc['i_ref']
        i_act = timeDc.get('i_dc_out', None)
    elif isinstance(timeAc, dict) and 'i_ref' in timeAc:
        # Use DC reference if present, otherwise AC
        if isinstance(timeAc['i_ref'], dict):
            if 'i_dc_out' in timeAc['i_ref']:
                i_ref = timeAc['i_ref']['i_dc_out']
                i_act = timeDc.get('i_dc_out', None)
            else:
                i_ref = timeAc['i_ref'].get('A', None)
                i_act = timeAc.get('i_ac_pri', None)
        else:
            i_ref = timeAc['i_ref']
            i_act = timeAc.get('i_ac_pri', None)
    else:
        i_act = timeDc.get('i_dc_out', timeAc.get('i_ac_pri', None))

    if i_act is not None:
        t_i = _time_for(i_act)
        plt.figure(figsize=(12, 6))
        plt.title("Closed Loop DAB: Reference Tracking")
        if i_ref is not None:
            t_r = _time_for(i_ref)
            plt.plot(t_r[::down2], i_ref[::down2], label="Reference")
        plt.plot(t_i[::down2], i_act[::down2], label="Actual")
        if i_ref is not None:
            err = i_act - i_ref
            plt.plot(t_i[::down2], err[::down2], label="Error")
        plt.xlabel("Time (s)")
        plt.ylabel("Current (A)")
        plt.grid('on')
        plt.legend(loc='upper right')

    ###################################################################################################################
    # DC-Link Voltages
    ###################################################################################################################
    if 'v_dc_pri' in timeDc and 'v_dc_sec' in timeDc:
        t_v = _time_for(timeDc['v_dc_pri'])
        plt.figure(figsize=(12, 6))
        plt.title("Closed Loop DAB: DC-Link Voltages")
        plt.plot(t_v[::down2], timeDc['v_dc_pri'][::down2], label="Primary $v_{dc}$")
        plt.plot(t_v[::down2], timeDc['v_dc_sec'][::down2], label="Secondary $v_{dc}$")
        plt.xlabel("Time (s)")
        plt.ylabel("Voltage (V)")
        plt.grid('on')
        plt.legend(loc='upper right')

    ###################################################################################################################
    # Switching States (short window)
    ###################################################################################################################
    def _get_col(obj, key):
        try:
            return obj[key]
        except Exception:
            return None

    sA = _get_col(timeSw, 'sA')
    sC = _get_col(timeSw, 'sC')
    t_s_col = _get_col(timeSw, 't')
    if sA is not None and sC is not None and t_s_col is not None:
        t_s = t_s_col.values if hasattr(t_s_col, 'values') else t_s_col
        plt.figure(figsize=(12, 4))
        plt.title("Closed Loop DAB: Switching States")
        plt.plot(t_s[::down2], sA[::down2], label="$s_A$")
        plt.plot(t_s[::down2], sC[::down2], label="$s_C$")
        plt.xlabel("Time (s)")
        plt.ylabel("State (p.u.)")
        plt.grid('on')
        plt.legend(loc='upper right')

    plt.show()

    ###################################################################################################################
    # MSG OUT
    ###################################################################################################################
    print("END: Plotting Closed Loop DAB Waveforms")
