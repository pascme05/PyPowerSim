#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         classB2
# Date:         27.04.2024
# Author:       Dr. Pascal A. Schirmer
# Version:      V.1.0
# Copyright:    Pascal Schirmer
#######################################################################################################################
#######################################################################################################################

#######################################################################################################################
# Function Description
#######################################################################################################################
"""
This class initialises an object of the B2 (half-bridge) circuit.
Inputs:     1) fel:     electrical frequency at the output (Hz)
            2) fs:      switching frequency of the converter cell (Hz)
            3) fsim:    simulation frequency of the converter cell (Hz)
            4) td:      dead-time of the switching devices (sec)
            5) tmin:    minimum on-time of a pulse (sec)
            6) cyc:     number of cycles till convergence
            7) W:       number of points evaluated for distortion analysis
            8) Mi:      modulation index (p.u.)
            9) Vdc:     dc link voltage at the input of the converter cell (V)
            10) Tc_st:  case temperature for stationary analysis
            11) Tj_st:  core temperature for stationary analysis
            12) Tc_tr:  case temperature for transient analysis
            13) Tj_tr:  core temperature for transient analysis
"""

#######################################################################################################################
# Import libs
#######################################################################################################################
# ==============================================================================
# Internal
# ==============================================================================
from src.general.helpFnc import cbInter, con2dis, deadTime
from src.pwm.oppPWM import oppPWM
from src.pwm.genWaveform import genWave
from src.cont.conHys import conHys

# ==============================================================================
# External
# ==============================================================================
import numpy as np
import pandas as pd
from scipy import signal
import cmath


#######################################################################################################################
# Class
#######################################################################################################################
class classB2:
    ###################################################################################################################
    # Constructor
    ###################################################################################################################
    def __init__(self, fel, fs, fsim, td, tmin, cyc, W, Mi, Vdc, Tc_st, Tj_st, Tc_tr, Tj_tr):
        self.fel = fel
        self.fs = fs
        self.fsim = fsim
        self.td = td
        self.tmin = tmin
        self.K = int(cyc)
        self.W = int(W)
        self.Mi = Mi
        self.Vdc = Vdc
        self.Tc_st = Tc_st
        self.Tj_st = Tj_st
        self.Tc_tr = Tc_tr
        self.Tj_tr = Tj_tr
        self.Nsim = int(np.ceil(fsim/fel))
        self.Npwm = int(np.ceil(fs/fel))
        self.q = int(fs / fel)
        self.N = int(fsim / fel)
        self.Ts = 1 / fs
        self.Tel = 1 / fel
        self.Mi_max = 1
        self.id1 = ['A']
        self.id2 = ['S1', 'S2']
        self.id3 = ['A', 'A']
        self.id4 = ['i_a', 'i_a']
        self.id5 = ['HS', 'LS']
        self.id6 = ['T1', 'T2']
        self.id7 = ['D1', 'D2']
        self.id8 = ['C1', 'C2']
        self.id9 = [1, 1]
        self.name = 'B2'

    ###################################################################################################################
    # Init Data
    ###################################################################################################################
    def initData(self):
        # ==============================================================================
        # Description
        # ==============================================================================
        """
        This function initialises the data structure of the B2 bridge. This method is
        supposed to be static.

        Input:

        Output:
        1) data:        initialsied data structure
        """

        # ==============================================================================
        # Initialisation
        # ==============================================================================
        data = {'elec': {}, 'loss': {}, 'ther': {}}

        # ==============================================================================
        # Switches
        # ==============================================================================
        # ------------------------------------------
        # Electric
        # ------------------------------------------
        data['elec']['sw'] = {}
        data['elec']['sw']['S1'] = pd.DataFrame(columns=['i_T', 'v_T', 'i_D', 'v_D'])
        data['elec']['sw']['S2'] = pd.DataFrame(columns=['i_T', 'v_T', 'i_D', 'v_D'])

        # ------------------------------------------
        # Losses
        # ------------------------------------------
        data['loss']['sw'] = {}
        data['loss']['sw']['S1'] = pd.DataFrame(columns=['p_T_c', 'p_T_s', 'p_T', 'p_D_c', 'p_D_s', 'p_D', 'p_L'])
        data['loss']['sw']['S2'] = pd.DataFrame(columns=['p_T_c', 'p_T_s', 'p_T', 'p_D_c', 'p_D_s', 'p_D', 'p_L'])

        # ------------------------------------------
        # Thermal
        # ------------------------------------------
        data['ther']['sw'] = pd.DataFrame(columns=['T1', 'T2', 'D1', 'D2', 'C1', 'C2'])

        # ==============================================================================
        # Capacitor
        # ==============================================================================
        # ------------------------------------------
        # Electric
        # ------------------------------------------
        data['elec']['cap'] = {}
        data['elec']['cap']['C1'] = pd.DataFrame(columns=['i_c', 'v_c'])

        # ------------------------------------------
        # Losses
        # ------------------------------------------
        data['loss']['cap'] = {}
        data['loss']['cap']['C1'] = pd.DataFrame(columns=['p_L'])

        # ------------------------------------------
        # Thermal
        # ------------------------------------------
        data['ther']['cap'] = pd.DataFrame(columns=['C1'])

        # ==============================================================================
        # Return
        # ==============================================================================
        return data

    ###################################################################################################################
    # Init Output
    ###################################################################################################################
    def initOut(self):
        # ==============================================================================
        # Description
        # ==============================================================================
        """
        This function initialises the output files of the B2 bridge

        Input:

        Output:
        1) timeSw:      time domain switching function
        2) timeElec:    time domain electrical
        3) timeLoss:    time domain losses
        4) timeTher:    time domain thermal
        5) freqSw:      frequency domain switching function
        6) freqDc:      frequency domain dc
        7) freqAc:      frequency domain ac
        8) distAc:      distortion ac
        9) distDc:      distortion dc
        """

        # ==============================================================================
        # Initialisation
        # ==============================================================================
        distAc = {}
        distDc = {}
        timeElec = {'sw': {}, 'cap': {}}
        timeLoss = {'sw': {}, 'cap': {}}
        timeTher = {'sw': {}, 'cap': {}}

        # ==============================================================================
        # Calculation
        # ==============================================================================
        timeElec['cap']['C1'] = pd.DataFrame(columns=['i_c', 'v_c'])
        timeSw = pd.DataFrame(columns=['t', 'v_a_ref', 'e_a', 'xAs', 'xAsh', 'sA', 'c'])
        freqSw = pd.DataFrame(columns=['S', 'Xs'])
        freqAc = pd.DataFrame(columns=['I_a', 'V_a'])
        freqDc = pd.DataFrame(columns=['I_dc', 'I_d_p', 'I_d_m', 'V_dc'])
        distAc['num'] = pd.DataFrame(data=np.ones((self.W, 6)),
                                     columns=['V_a_eff', 'V_a_v1_eff', 'V_a_thd', 'I_a_eff', 'I_a_v1_eff', 'I_a_thd'])
        distAc['ana'] = pd.DataFrame(data=np.ones((self.W, 6)),
                                     columns=['V_a_eff', 'V_a_v1_eff', 'V_a_thd', 'I_a_eff', 'I_a_v1_eff', 'I_a_thd'])
        distDc['num'] = pd.DataFrame(data=np.ones((self.W, 6)),
                                     columns=['V_dc_eff', 'V_dc_v1_eff', 'V_dc_thd', 'I_dc_eff', 'I_dc_v1_eff', 'I_dc_thd'])
        distDc['ana'] = pd.DataFrame(data=np.ones((self.W, 6)),
                                     columns=['V_dc_eff', 'V_dc_v1_eff', 'V_dc_thd', 'I_dc_eff', 'I_dc_v1_eff', 'I_dc_thd'])

        # ==============================================================================
        # Return
        # ==============================================================================
        return [timeSw, timeElec, timeLoss, timeTher, freqSw, freqDc, freqAc, distAc, distDc]

    ###################################################################################################################
    # Init Output
    ###################################################################################################################
    def out(self, timeElec, timeLoss, timeTher, timeAc, timeDc, freqSw, freqAc, freqDc, distAc, distDc, t_ref, v_ref, e_ref, s, c, xs, xsh, x, xN0, M_i, t0, t1, Nel):
        # ==============================================================================
        # Description
        # ==============================================================================
        """
        This function summarizes the output structure of the B2 bridge.

        Input:
        1) timeElec:    temporal electrical signals
        2) timeLoss:    temporal loss signals
        3) timeTher:    temporal thermal signals
        4) timeAc:      ac time domain signals
        5) timeDc:      dc time domain signals
        6) freqSw:      frequency domain switching functions
        7) freqAc:      freq domain switching functions ac
        8) freqDc:      freq domain switching functions dc
        8) distAc:      distortion domain ac signals
        9) distDc:      distortion domain dc signals
        10) t_ref:      reference time vector (sec)
        11) v_ref:      reference voltage vector (V)
        12) e_ref:      reference back emf vector (V)
        13) s:          switching function time domain
        14) c:          carrier waveform time domain
        15) xs:         sampled reference signal
        16) xsh:        sampled reference signal including zero order hold
        17) x:          reference signal
        18) xN0:        reference zero sequence
        19) M_i:        vector of modulation indices
        20) t0:         starting sample
        21) t1:         ending sample
        22) Nel:        number of electrical cycles

        Output:
        1) time:        time domain output signals
        2) freq:        frequency domain output signals
        3) sweep:       sweeping domain output signals

        """

        # ==============================================================================
        # Initialisation
        # ==============================================================================
        time = {}
        freq = {}
        sweep = {}
        [timeSw, _, _, _, _, _, _, _, _] = self.initOut()

        # ==============================================================================
        # Calculation
        # ==============================================================================
        timeSw['t'] = t_ref[0:(t1 - t0)]
        timeSw['v_ref'] = v_ref['A'][t0:t1]
        timeSw['e'] = e_ref['A'][t0:t1]
        timeSw['sA'] = s['A'][t0:t1]
        timeSw['c'] = c['A'][t0:t1]
        timeSw['xAs'] = xs['A'][t0:t1]
        timeSw['xAsh'] = xsh['A'][t0:t1]
        timeSw['xA'] = x['A'][t0:t1]
        timeSw['n0'] = xN0['A'][t0:t1]

        # ==============================================================================
        # Combine
        # ==============================================================================
        # ------------------------------------------
        # Time
        # ------------------------------------------
        # Transient time
        if not timeLoss:
            time['t'] = timeSw['t']
        else:
            time['t'] = np.linspace(0, self.Tel * Nel, int(len(timeLoss['sw']['S1']['p_T'])))

        # Variables
        time['Sw'] = timeSw
        time['Ac'] = timeAc
        time['Dc'] = timeDc
        time['Elec'] = timeElec
        time['Loss'] = timeLoss
        time['Ther'] = timeTher

        # ------------------------------------------
        # Frequency
        # ------------------------------------------
        freq['Sw'] = freqSw
        freq['Ac'] = freqAc
        freq['Dc'] = freqDc

        # ------------------------------------------
        # Sweep
        # ------------------------------------------
        sweep['Ac'] = distAc
        sweep['Dc'] = distDc
        sweep['Mi'] = M_i

        # ==============================================================================
        # Return
        # ==============================================================================
        return [time, freq, sweep]

    ###################################################################################################################
    # Dead Time and Minimum Pulse Width
    ###################################################################################################################
    def calcDead(self, s, t_ref, Mi):
        # ==============================================================================
        # Description
        # ==============================================================================
        """
        This function calculates the dead-time and the minimum pulse-width on a switching
        sequence.

        Input:
        1) s:       Input switching sequence
        2) t_ref:   Reference time (sec)
        3) Mi:      Modulation index (0 ... 4/pi)

        Output:
        1) s:      updated switching sequence
        """

        # ==============================================================================
        # Initialisation
        # ==============================================================================
        tmin = int(self.tmin / (t_ref[1] - t_ref[0]))
        td = int(self.td / (t_ref[1] - t_ref[0]))

        # ==============================================================================
        # Calculation
        # ==============================================================================
        # ------------------------------------------
        # Minimum Pulse Width
        # ------------------------------------------
        if self.tmin > 0:
            hold = tmin
            for i in range(0, len(s)):
                if hold >= tmin:
                    if Mi != 0:
                        s[i] = s[i]
                    hold = 0
                else:
                    s[i] = s[i - 1]
                    hold = hold + 1

        # ------------------------------------------
        # Dead-time
        # ------------------------------------------
        if self.td > 0:
            s = deadTime(s, td)

        # ==============================================================================
        # Return
        # ==============================================================================
        return s

    ###################################################################################################################
    # Fundamental Frequency Control
    ###################################################################################################################
    def calcSeqFF(self, v_ref, t_ref, Mi):
        # ==============================================================================
        # Description
        # ==============================================================================
        """
        This function calculates the switching time assuming fundamental frequency
        control.

        Input:
        1) v_ref:   Reference voltage (V)
        2) t_ref:   Reference time (sec)
        3) Mi:      Modulation index (0 ... 4/pi)

        Output:
        1) xs:      sampled reference signal
        2) xsh:     sampled reference signal including zero order hold
        3) s:       switching instances
        4) c:       carrier signal
        5) x:       reference signal
        6) xN0:     zero sequence reference signal
        """

        # ==============================================================================
        # Initialisation
        # ==============================================================================
        c = np.zeros(np.size(t_ref))

        # ==============================================================================
        # Calculation
        # ==============================================================================
        x = Mi * v_ref / np.max(v_ref)
        xN0 = np.zeros(np.size(x))
        s = signal.square(2 * np.pi * self.fel * t_ref, duty=0.5)
        xs = x
        xsh = x

        # ==============================================================================
        # Return
        # ==============================================================================
        return [xs, xsh, s, c, x, xN0]

    ###################################################################################################################
    # Carrier based PWM
    ###################################################################################################################
    def calcSeqCB(self, v_ref, t_ref, Mi, setup):
        # ==============================================================================
        # Description
        # ==============================================================================
        """
        This function calculates the switching time assuming carrier based PWM.

        Input:
        1) v_ref:   Reference voltage (V)
        2) t_ref:   Reference time (sec)
        3) Mi:      Modulation index (0 ... 4/pi)

        Output:
        1) xs:      sampled reference signal
        2) xsh:     sampled reference signal including zero order hold
        3) s:       switching instances
        4) c:       carrier signal
        5) x:       reference signal
        6) xN0:     zero sequence reference signal
        """

        # ==============================================================================
        # Initialisation
        # ==============================================================================
        tmin = int(self.tmin / (t_ref[1] - t_ref[0]))

        # ==============================================================================
        # Calculation
        # ==============================================================================
        # ------------------------------------------
        # Reference
        # ------------------------------------------
        x = Mi * v_ref / np.max(v_ref)
        xN0 = np.zeros(np.size(x))

        # ------------------------------------------
        # Carrier
        # ------------------------------------------
        if setup['Par']['PWM']['tri'] == "RE":
            c = signal.sawtooth(2 * np.pi * self.fs * t_ref, 1) * (-1)
        elif setup['Par']['PWM']['tri'] == "FE":
            c = signal.sawtooth(2 * np.pi * self.fs * (t_ref - 0.5 / self.fs), 0) * (-1)
        elif setup['Par']['PWM']['tri'] == "AM":
            c = signal.sawtooth(2 * np.pi * self.fs * t_ref, 1 / 3) * (-1)
        else:
            c = signal.sawtooth(2 * np.pi * self.fs * t_ref, 0.5) * (-1)
        c = (2 * (c - min(c)) / (max(c) - min(c))) - 1

        # ------------------------------------------
        # Sampling
        # ------------------------------------------
        if setup['Par']['PWM']['samp'] == "RS":
            if setup['Par']['PWM']['upd'] == "SE":
                xs = con2dis(x, t_ref, self.Ts)
                xsh = np.roll(x, int(len(xs) * self.fel / self.fs))
            else:
                xs = con2dis(x, t_ref, self.Ts / 2)
                xsh = np.roll(x, int(len(xs) * self.fel / self.fs / 2))
        else:
            xs = x
            xsh = x

        # ------------------------------------------
        # Intersections
        # ------------------------------------------
        s = cbInter(xs, c, Mi, tmin)
        s = self.calcDead(s, t_ref, Mi)

        # ==============================================================================
        # Return
        # ==============================================================================
        return [xs, xsh, s, c, x, xN0]

    ###################################################################################################################
    # Optimal Pulse Patterns
    ###################################################################################################################
    def calcSeqOPP(self, v_ref, t_ref, Mi, setup):
        # ==============================================================================
        # Description
        # ==============================================================================
        """
        This function calculates the switching time assuming optimal pulse patterns.

        Input:
        1) v_ref:   Reference voltage (V)
        2) t_ref:   Reference time (sec)
        3) Mi:      Modulation index (0 ... 4/pi)

        Output:
        1) xs:      sampled reference signal
        2) xsh:     sampled reference signal including zero order hold
        3) s:       switching instances
        4) c:       carrier signal
        5) x:       reference signal
        6) xN0:     zero sequence reference signal
        """

        # ==============================================================================
        # Initialisation
        # ==============================================================================
        kmax = 10 * self.q
        ang = np.linspace(0, np.pi * 2 * self.N, np.size(t_ref))
        s = (-1) * np.ones(np.size(t_ref))
        c = np.zeros(np.size(t_ref))
        ang_total = []
        val_total = []

        # ==============================================================================
        # Calculation
        # ==============================================================================
        # ------------------------------------------
        # Reference
        # ------------------------------------------
        x = Mi * v_ref / np.max(v_ref)
        xN0 = np.zeros(np.size(x))

        # ------------------------------------------
        # Optimal Angles
        # ------------------------------------------
        # Fundamental angle (0 ... 2pi)
        [ang_fun, val_fun, _] = oppPWM(kmax, self.q, Mi / 4 * np.pi, 4, setup)

        # Complete angle
        for i in range(0, self.N):
            if i == 0:
                ang_total = ang_fun
                val_total = val_fun
            else:
                ang_total = np.concatenate((ang_total, ang_fun + 2 * np.pi * i), axis=0)
                val_total = np.concatenate((val_total, val_fun), axis=0)

        # ------------------------------------------
        # Switching times
        # ------------------------------------------
        # Switching Edges
        for i in range(0, len(ang_total)):
            idx = np.argmin(abs(ang - ang_total[i]))
            c[idx] = (-1) ** i

        # Switching Function
        for i in range(1, len(s)):
            if c[i] == 0:
                s[i] = s[i - 1]
            else:
                s[i] = c[i]

        # Direction
        for i in range(1, len(s)):
            if v_ref[i] >= 0:
                s[i] = s[i] * (-1)
            else:
                s[i] = s[i] * (+1)

        # Dead time
        s = self.calcDead(s, t_ref, Mi)

        # ------------------------------------------
        # Sampling
        # ------------------------------------------
        if setup['Par']['PWM']['samp'] == "RS":
            if setup['Par']['PWM']['upd'] == "SE":
                xs = con2dis(x, t_ref, self.Ts)
                xsh = np.roll(x, int(len(xs) * self.fel / self.fs))
            else:
                xs = con2dis(x, t_ref, self.Ts / 2)
                xsh = np.roll(x, int(len(xs) * self.fel / self.fs / 2))
        else:
            xs = x
            xsh = x

        # ==============================================================================
        # Return
        # ==============================================================================
        return [xs, xsh, s, c, x, xN0]

    ###################################################################################################################
    # Carrier based PWM
    ###################################################################################################################
    def calcPWM(self, v_ref, t_ref, Mi, setup):
        # ==============================================================================
        # Description
        # ==============================================================================
        """
        This function calculates the switching times assuming any provided PWM method. Valid
        PWM methods are carrier-based, fundamental frequency, and optimal pulse patterns.

        Input:
        1) v_ref:   Reference voltage (V)
        2) t_ref:   Reference time (sec)
        3) Mi:      Modulation index (0 ... 4/pi)

        Output:
        1) xs:      sampled reference signal
        2) xsh:     sampled reference signal including zero order hold
        3) s:       switching instances
        4) c:       carrier signal
        5) x:       reference signal
        6) xN0:     zero sequence reference signal
        """

        # ==============================================================================
        # Init
        # ==============================================================================
        x_out = {}
        xs_out = {}
        xsh_out = {}
        c_out = {}
        s_out = {}
        xN0_out = {}

        # ==============================================================================
        # Calculation
        # ==============================================================================
        if setup['Par']['PWM']['type'] == "FF":
            [xs, xsh, s, c, x, xN0] = self.calcSeqFF(v_ref['A'], t_ref, Mi)
        elif setup['Par']['PWM']['type'] == "CB":
            [xs, xsh, s, c, x, xN0] = self.calcSeqCB(v_ref['A'], t_ref, Mi, setup)
        elif setup['Par']['PWM']['type'] == "OPP":
            [xs, xsh, s, c, x, xN0] = self.calcSeqOPP(v_ref['A'], t_ref, Mi, setup)
        else:
            [xs, xsh, s, c, x, xN0] = self.calcSeqCB(v_ref['A'], t_ref, Mi, setup)

        # ==============================================================================
        # Post-Processing
        # ==============================================================================
        xs_out['A'] = xs
        xsh_out['A'] = xsh
        s_out['A'] = s
        c_out['A'] = c
        x_out['A'] = x
        xN0_out['A'] = xN0

        # ==============================================================================
        # Return
        # ==============================================================================
        return [xs_out, xsh_out, s_out, c_out, x_out, xN0_out]

    ###################################################################################################################
    # Closed Loop Control
    ###################################################################################################################
    def calcCON(self, i_ref, i_act, s_act, setup):
        # ==============================================================================
        # Description
        # ==============================================================================
        """
        This function calculates the switching times assuming closed loop control. It
        compares the reference current and the actual current and updates the switching
        state accordingly using an arbitrary control algorithm.

        Input:
        1) i_ref:   Reference current (A)
        2) i_act:   Actual current (A)
        3) s_act:   Actual switching states
        4) setup:   variable including all parameters

        Output:
        1) s:       switching instances (sec)
        2) Mi:      updated modulation index
        3) err:     error between reference and actual signal
        """

        # ==============================================================================
        # Init
        # ==============================================================================
        s_out = {}
        tol = np.max(abs(i_ref)) * setup['Par']['Cont']['hys'] / 100

        # ==============================================================================
        # Calculation
        # ==============================================================================
        if setup['Par']['Cont']['type'] == "HY":
            [s, err] = conHys(i_act, i_ref, s_act['A'], tol)
        elif setup['Par']['Cont']['type'] == "PI":
            [s, err] = conHys(i_act, i_ref, s_act['A'], tol)
        else:
            [s, err] = conHys(i_act, i_ref, s_act['A'], tol)

        # ==============================================================================
        # Post-Processing
        # ==============================================================================
        s_out['A'] = s
        Mi = 4 / np.pi * abs(np.mean(s))

        # ==============================================================================
        # Return
        # ==============================================================================
        return [s_out, Mi, err]

    ###################################################################################################################
    # Temporal Output
    ###################################################################################################################
    def calcTime(self, s, e, t, Mi, mdl, t0, t1, init, avg, setup):
        # ==============================================================================
        # Description
        # ==============================================================================
        """
        This function calculates the temporal outputs of the B2 bridge (AC and DC side).

        Input:
        1) s:       switching function
        2) e:       induced voltage
        3) t:       reference time (sec)
        4) Mi:      modulation index (0 ... 4/pi)
        5) mdl:     transfer functions
        6) t0:      start time (sample)
        7) t1:      end time (sample)
        8) setup:   all setup variables
        9) init:    initial conditions for the lsim solver
        10) avg:    if 1) mean is subtracted from the signal

        Output:
        1) outAc:   outputs ac side
        2) outDc:   outputs dc side
        3) outIn:   initial conditions for the next period
        """

        # ==============================================================================
        # Initialisation
        # ==============================================================================
        # ------------------------------------------
        # Outputs
        # ------------------------------------------
        outAc = {}
        outDc = {}
        outIn = {}

        # ------------------------------------------
        # Initial Conditions
        # ------------------------------------------
        if not init:
            init = {'inp': 0, 'out': 0, 'dc': 0, 'load': 0}

        # ==============================================================================
        # Calculation
        # ==============================================================================
        # ------------------------------------------
        # AC Side
        # ------------------------------------------
        # Inverter Output
        v_a0 = 0.5 * s['A'] * self.Vdc

        # Filter Output
        if setup['Top']['outFilter'] == 0:
            v_L = v_a0
        else:
            _, v_L, _, = signal.lsim(mdl['SS']['Out'], v_a0, t, X0=init['out'])

        # Load
        v_a = v_L - Mi * e['A']

        # Current
        if setup['Top']['wave'] == "con":
            _, i_a, _, = signal.lsim(mdl['SS']['Load'], v_a, t, X0=init['load'])
            i_a = i_a[t0:t1]
        else:
            if avg == 1:
                _, i_a, _, = signal.lsim(mdl['SS']['Load'], (v_a - np.mean(v_a)), t, X0=init['load'])
                i_a = i_a[t0:t1]
                i_a = i_a - np.mean(i_a)
            else:
                _, i_a, _, = signal.lsim(mdl['SS']['Load'], v_a, t, X0=init['load'])
                i_a = i_a[t0:t1]

        # ------------------------------------------
        # DC Side
        # ------------------------------------------
        # Inverter Input
        i_d_p = i_a * (1 + s['A'][t0:t1]) / 2
        i_d_m = i_a * (1 - s['A'][t0:t1]) / 2
        i_dc = i_d_p

        # DC-Link
        i_c = np.mean(i_d_p) - i_d_p
        _, v_dc, _, = signal.lsim(mdl['SS']['DC'], i_c, t[t0:t1], X0=init['dc'])

        # Filter Input
        if setup['Top']['inpFilter'] == 0:
            v_in = v_dc
        else:
            _, v_in, _, = signal.lsim(mdl['SS']['Inp'], (v_dc - self.Vdc), t[t0:t1], X0=init['inp'])
            v_in = v_in + self.Vdc

        # ==============================================================================
        # Post-Processing
        # ==============================================================================
        # ------------------------------------------
        # AC Side
        # ------------------------------------------
        outAc['v_a0'] = v_a0[t0:t1]
        outAc['v_L'] = v_L[t0:t1]
        outAc['v_a_out'] = v_L[t0:t1]
        outAc['v_a'] = v_a[t0:t1]
        outAc['i_a'] = i_a

        # ------------------------------------------
        # DC Side
        # ------------------------------------------
        outDc['v_in'] = v_in
        outDc['v_dc'] = v_dc
        outDc['i_dc'] = i_dc
        outDc['i_c'] = i_c
        outDc['i_d_m'] = i_d_m
        outDc['i_d_p'] = i_d_p

        # ------------------------------------------
        # Init Conditions
        # ------------------------------------------
        outIn['inp'] = v_in[-1]
        outIn['out'] = v_L[-1]
        outIn['dc'] = v_dc[-1]
        outIn['load'] = i_a[-1]

        # ==============================================================================
        # Return
        # ==============================================================================
        return [outAc, outDc, outIn]

    ###################################################################################################################
    # Analytical distortion
    ###################################################################################################################
    def calcDist(self, _1, _2, Mi, L, Z, _3):
        # ==============================================================================
        # Description
        # ==============================================================================
        """
        This function calculates the analytical distortion of the B2 bridge.

        Input:
        1) i_a:     phase current (A)
        2) v_a:     phase voltage (V)
        3) Mi:      modulation index (0 ... 4/pi)
        4) L:       load inductance
        5) Z:       load impedance
        6) setup:   variable including all setup parameters

        Output:
        1) outAc:   outputs distortion ac side
        2) outDc:   outputs distortion dc side
        """

        # ==============================================================================
        # Initialisation
        # ==============================================================================
        outAc = {}
        outDc = {}

        # ==============================================================================
        # Calculation
        # ==============================================================================
        # ------------------------------------------
        # AC Side
        # ------------------------------------------
        V_a0_eff = self.Vdc / 2
        V_a0_v1_eff = 1 / np.sqrt(2) * self.Vdc / 2 * Mi
        V_a0_thd = self.Vdc / 2 * np.sqrt(1 - Mi ** 2 / 2)
        I_a_eff = V_a0_eff / Z
        I_a_v1_eff = V_a0_v1_eff / Z
        I_a_thd = 1 / np.sqrt(48) * self.Vdc / 2 * self.Ts / L * np.sqrt(3 / 8 * Mi ** 4 - Mi ** 2 + 1)

        # ------------------------------------------
        # DC Side
        # ------------------------------------------
        V_dc_eff = Mi * 4 / np.pi * self.Vdc / 2

        # ==============================================================================
        # Post-Processing
        # ==============================================================================
        # ------------------------------------------
        # AC Side
        # ------------------------------------------
        outAc['V_a_eff'] = V_a0_eff
        outAc['V_a_v1_eff'] = V_a0_v1_eff
        outAc['V_a_thd'] = V_a0_thd
        outAc['I_a_eff'] = I_a_eff
        outAc['I_a_v1_eff'] = I_a_v1_eff
        outAc['I_a_thd'] = I_a_thd

        # ------------------------------------------
        # DC Side
        # ------------------------------------------
        outDc['V_dc_eff'] = V_dc_eff
        outDc['V_dc_v1_eff'] = np.nan
        outDc['V_dc_thd'] = np.nan
        outDc['I_dc_eff'] = np.nan
        outDc['I_dc_v1_eff'] = np.nan
        outDc['I_dc_thd'] = np.nan

        # ==============================================================================
        # Return
        # ==============================================================================
        return [outAc, outDc]

    ###################################################################################################################
    # Calculations frequency domain
    ###################################################################################################################
    def calcRef(self, E, phiE, phiV, t, setup):
        # ==============================================================================
        # Description
        # ==============================================================================
        """
        This function calculates the reference voltage, the reference current and back EMF
        functions based on the B2 topology and the given parameters.

        Input:
        1) E:       amplitude of the back emf (V)
        2) phiE:    angle of the back emf (rad)
        3) v_a:     load angle of the output (rad)
        4) t:       given time vector (sec)
        5) setup:   file including all setup parameters

        Output:
        1) v_ref:   reference voltage for given load scenario (V)
        2) e_ref:   reference back emf for given load scenario (V)
        3) i_ref:   reference current for given load scenario (A)
        """

        # ==============================================================================
        # Initialisation
        # ==============================================================================
        v_ref = {}
        i_ref = {}
        e_ref = {}
        Io = cmath.polar(setup['Dat']['stat']['Io'])[0]

        # ==============================================================================
        # Calculation
        # ==============================================================================
        # ------------------------------------------
        # Time
        # ------------------------------------------
        try:
            if not t:
                t = np.linspace(0, self.K / self.fel, self.K * self.N + 1)
        except:
            test = 1

        # ------------------------------------------
        # Reference
        # ------------------------------------------
        v_ref['A'] = (self.Vdc / 2) * self.Mi * genWave(t, self.fel, phiV, setup)
        e_ref['A'] = E * genWave(t, self.fel, phiE, setup)
        i_ref['A'] = Io * np.sqrt(2) * genWave(t, self.fel, setup['Dat']['stat']['PhiVI'], setup)

        # ==============================================================================
        # Return
        # ==============================================================================
        return [v_ref, e_ref, i_ref]

#######################################################################################################################
# References
#######################################################################################################################
