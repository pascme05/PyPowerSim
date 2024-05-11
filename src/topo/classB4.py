#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         classB4
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
This class initialises an object of the B4 (full-bridge) circuit.
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

# ==============================================================================
# External
# ==============================================================================
import numpy as np
import pandas as pd
from scipy import signal


#######################################################################################################################
# Class
#######################################################################################################################
class classB4:
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
        self.id1 = ['A', 'B']
        self.id2 = ['S1', 'S2', 'S3', 'S4']
        self.id3 = ['A', 'A', 'B', 'B']
        self.id4 = ['i_a', 'i_a', 'i_a', 'i_a']
        self.id5 = ['HS', 'LS', 'HS', 'LS']
        self.id6 = ['T1', 'T2', 'T3', 'T4']
        self.id7 = ['D1', 'D2', 'D3', 'D4']
        self.id8 = ['C1', 'C2', 'C3', 'C4']
        self.id9 = [+1, +1, -1, -1]
        self.name = 'B4'

    ###################################################################################################################
    # Init Data
    ###################################################################################################################
    def initData(self):
        # ==============================================================================
        # Description
        # ==============================================================================
        """
        This function initialises the data structure of the B4 bridge. This method is
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
        data['elec']['sw']['S3'] = pd.DataFrame(columns=['i_T', 'v_T', 'i_D', 'v_D'])
        data['elec']['sw']['S4'] = pd.DataFrame(columns=['i_T', 'v_T', 'i_D', 'v_D'])

        # ------------------------------------------
        # Losses
        # ------------------------------------------
        data['loss']['sw'] = {}
        data['loss']['sw']['S1'] = pd.DataFrame(columns=['p_T_c', 'p_T_s', 'p_T', 'p_D_c', 'p_D_s', 'p_D', 'p_L'])
        data['loss']['sw']['S2'] = pd.DataFrame(columns=['p_T_c', 'p_T_s', 'p_T', 'p_D_c', 'p_D_s', 'p_D', 'p_L'])
        data['loss']['sw']['S3'] = pd.DataFrame(columns=['p_T_c', 'p_T_s', 'p_T', 'p_D_c', 'p_D_s', 'p_D', 'p_L'])
        data['loss']['sw']['S4'] = pd.DataFrame(columns=['p_T_c', 'p_T_s', 'p_T', 'p_D_c', 'p_D_s', 'p_D', 'p_L'])

        # ------------------------------------------
        # Thermal
        # ------------------------------------------
        data['ther']['sw'] = pd.DataFrame(columns=['T1', 'T2', 'T3', 'T4',
                                                   'D1', 'D2', 'D3', 'D4',
                                                   'C1', 'C2', 'C3', 'C4'])

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
        This function initialises the output files of the B4 bridge

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
        timeSw = pd.DataFrame(columns=['t', 'v_a_ref', 'v_b_ref', 'e', 'xAs', 'xBs', 'xAsh', 'xBsh', 'sA', 'sB', 'cA', 'cB'])
        freqSw = pd.DataFrame(columns=['Sa', 'Sb', 'Xas', 'Xbs'])
        freqAc = pd.DataFrame(columns=['I_ab', 'V_ab'])
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
        This function summarizes the output structure of the B4 bridge.

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
        timeSw['v_a_ref'] = v_ref['A'][t0:t1]
        timeSw['v_b_ref'] = v_ref['B'][t0:t1]
        timeSw['e_a'] = e_ref['A'][t0:t1]
        timeSw['sA'] = s['A'][t0:t1]
        timeSw['sB'] = s['B'][t0:t1]
        timeSw['c'] = c['A'][t0:t1]
        timeSw['cA'] = c['A'][t0:t1]
        timeSw['cB'] = c['B'][t0:t1]
        timeSw['xAs'] = xs['A'][t0:t1]
        timeSw['xBs'] = xs['B'][t0:t1]
        timeSw['xAsh'] = xsh['A'][t0:t1]
        timeSw['xBsh'] = xsh['B'][t0:t1]
        timeSw['xA'] = x['A'][t0:t1]
        timeSw['xB'] = x['B'][t0:t1]
        timeSw['n0'] = xN0['A'][t0:t1]
        timeSw['n0A'] = xN0['A'][t0:t1]
        timeSw['n0B'] = xN0['B'][t0:t1]

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
        c = {}
        s = {}
        x = {}
        alpha = np.arccos(Mi * np.pi / 4)

        # ==============================================================================
        # Calculation
        # ==============================================================================
        # ------------------------------------------
        # Reference
        # ------------------------------------------
        x['A'] = Mi * v_ref['A'] / np.max(np.abs(v_ref['A']))
        x['B'] = Mi * v_ref['B'] / np.max(np.abs(v_ref['B']))

        # ------------------------------------------
        # Carrier
        # ------------------------------------------
        c['A'] = np.zeros(np.size(t_ref))
        c['B'] = np.zeros(np.size(t_ref))

        # ------------------------------------------
        # Sequence
        # ------------------------------------------
        s['A'] = signal.square(2 * np.pi * self.fel * t_ref + np.pi / 2 - alpha, duty=0.5)
        s['B'] = signal.square(2 * np.pi * self.fel * t_ref + np.pi / 2 + alpha, duty=0.5)

        # ------------------------------------------
        # Output
        # ------------------------------------------
        xN0 = c
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
        # ------------------------------------------
        # Variables
        # ------------------------------------------
        x = {}
        c = {}
        s = {}
        xs = {}
        xsh = {}
        xN0 = {}

        # ------------------------------------------
        # Parameters
        # ------------------------------------------
        tmin = int(self.tmin / (t_ref[1] - t_ref[0]))

        # ==============================================================================
        # Calculation
        # ==============================================================================
        # ------------------------------------------
        # Reference
        # ------------------------------------------
        x['A'] = Mi * v_ref['A'] / np.max(np.abs(v_ref['A']))
        x['B'] = Mi * v_ref['B'] / np.max(np.abs(v_ref['B']))

        # ------------------------------------------
        # Carrier
        # ------------------------------------------
        # Interleaved
        if setup['Par']['PWM']['tri'] == "RE":
            c['A'] = signal.sawtooth(2 * np.pi * self.fs * t_ref, 1) * (-1)
            c['B'] = signal.sawtooth(2 * np.pi * self.fs * (t_ref - 0.5 / self.fs), 1)
        elif setup['Par']['PWM']['tri'] == "FE":
            c['A'] = signal.sawtooth(2 * np.pi * self.fs * (t_ref - 0.5 / self.fs), 0) * (-1)
            c['B'] = signal.sawtooth(2 * np.pi * self.fs * t_ref, 0)
        elif setup['Par']['PWM']['tri'] == "AM":
            c['A'] = signal.sawtooth(2 * np.pi * self.fs * t_ref, 1 / 3) * (-1)
            c['B'] = signal.sawtooth(2 * np.pi * self.fs * (t_ref - 0.5 / self.fs), 1 / 3)
        else:
            c['A'] = signal.sawtooth(2 * np.pi * self.fs * t_ref, 0.5) * (-1)
            c['B'] = signal.sawtooth(2 * np.pi * self.fs * (t_ref - 0.5 / self.fs), 0.5) * (-1)
        c['A'] = (2 * (c['A'] - min(c['A'])) / (max(c['A']) - min(c['A']))) - 1
        c['B'] = (2 * (c['B'] - min(c['B'])) / (max(c['B']) - min(c['B']))) - 1

        # Non-Interleaved
        if setup['Par']['PWM']['int'] == 0:
            c['B'] = c['A']

        # ------------------------------------------
        # Sampling
        # ------------------------------------------
        for i in range(0, len(self.id1)):
            if setup['Par']['PWM']['samp'] == "RS":
                if setup['Par']['PWM']['upd'] == "SE":
                    xs[self.id1[i]] = con2dis(x[self.id1[i]], t_ref, self.Ts)
                    xsh[self.id1[i]] = np.roll(x[self.id1[i]], int(len(xs[self.id1[i]]) * self.fel / self.fs))
                else:
                    xs[self.id1[i]] = con2dis(x[self.id1[i]], t_ref, self.Ts / 2)
                    xsh[self.id1[i]] = np.roll(x[self.id1[i]], int(len(xs[self.id1[i]]) * self.fel / self.fs / 2))
            else:
                xs[self.id1[i]] = x[self.id1[i]]
                xsh[self.id1[i]] = x[self.id1[i]]

        # ------------------------------------------
        # Intersections
        # ------------------------------------------
        s['A'] = cbInter(xs['A'], c['A'], Mi, tmin)
        s['B'] = cbInter(xs['B'], c['B'], Mi, tmin)

        # ------------------------------------------
        # Dead-Time
        # ------------------------------------------
        s['A'] = self.calcDead(s['A'], t_ref, Mi)
        s['B'] = self.calcDead(s['B'], t_ref, Mi)

        # ==============================================================================
        # Post-Processing
        # ==============================================================================
        xN0['A'] = np.zeros(np.size(c['A']))
        xN0['B'] = np.zeros(np.size(c['B']))

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
        # ------------------------------------------
        # Variables
        # ------------------------------------------
        x = {}
        c = {}
        s = {}
        xs = {}
        xsh = {}
        xN0 = {}

        # ------------------------------------------
        # Parameters
        # ------------------------------------------
        kmax = 10 * self.q
        ang = np.linspace(0, np.pi * 2 * self.N, np.size(t_ref))
        s['A'] = (-1) * np.ones(np.size(t_ref))
        c['A'] = np.zeros(np.size(t_ref))
        ss = np.zeros(np.size(t_ref))
        ang_total = []
        val_total = []

        # ==============================================================================
        # Calculation
        # ==============================================================================
        # ------------------------------------------
        # Reference
        # ------------------------------------------
        x['A'] = Mi * v_ref['A'] / np.max(np.abs(v_ref['A']))
        x['B'] = Mi * v_ref['B'] / np.max(np.abs(v_ref['B']))
        xN0['A'] = np.zero(np.size(x['A']))
        xN0['B'] = np.zero(np.size(x['B']))

        # ------------------------------------------
        # Optimal Angles
        # ------------------------------------------
        # Fundamental angle (0 ... 2pi)
        [ang_fun, val_fun, _] = oppPWM(kmax, self.q*2, Mi / 4 * np.pi, 4, setup)

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
            c['A'][idx] = (-1) ** i

        # Switching Function
        ang2 = ang
        for i in range(1, len(ss)):
            if ang2[i] < np.pi:
                if c['A'][i] == 0:
                    ss[i] = ss[i - 1]
                else:
                    if c['A'][i] > 0:
                        ss[i] = 1
                    else:
                        ss[i] = 0
            else:
                if c['A'][i] == 0:
                    ss[i] = ss[i - 1]
                else:
                    if c['A'][i] > 0:
                        ss[i] = -1
                    else:
                        ss[i] = 0

            if ang2[i] > 2 * np.pi:
                ang2 = ang2 - 2 * np.pi

        # Direction
        for i in range(1, len(ss)):
            if ss[i] >= 0:
                s['A'][i] = 1
            elif ss[i] < 0:
                s['A'][i] = -1

        # Interleaved
        s['B'] = np.roll(s['A'], -int(np.floor(180 / 360 / self.N * len(ss))))
        c['B'] = np.roll(c['A'], -int(np.floor(180 / 360 / self.N * len(ss))))

        # Non-Interleaved
        if setup['Par']['PWM']['int'] == 0:
            c['B'] = c['A']

        # Dead time
        s['A'] = self.calcDead(s['A'], t_ref, Mi)
        s['B'] = self.calcDead(s['B'], t_ref, Mi)

        # ------------------------------------------
        # Sampling
        # ------------------------------------------
        for i in range(0, len(self.id1)):
            if setup['Par']['PWM']['samp'] == "RS":
                if setup['Par']['PWM']['upd'] == "SE":
                    xs[self.id1[i]] = con2dis(x[self.id1[i]], t_ref, self.Ts)
                    xsh[self.id1[i]] = np.roll(x[self.id1[i]], int(len(xs[self.id1[i]]) * self.fel / self.fs))
                else:
                    xs[self.id1[i]] = con2dis(x[self.id1[i]], t_ref, self.Ts / 2)
                    xsh[self.id1[i]] = np.roll(x[self.id1[i]], int(len(xs[self.id1[i]]) * self.fel / self.fs / 2))
            else:
                xs[self.id1[i]] = x[self.id1[i]]
                xsh[self.id1[i]] = x[self.id1[i]]

        # ==============================================================================
        # Return
        # ==============================================================================
        return [xs, xsh, s, c, x, xN0]

    ###################################################################################################################
    # Calculate PWM
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
        # Calculation
        # ==============================================================================
        if setup['Par']['PWM']['type'] == "FF":
            [xs, xsh, s, c, x, xN0] = self.calcSeqFF(v_ref, t_ref, Mi)
        elif setup['Par']['PWM']['type'] == "CB":
            [xs, xsh, s, c, x, xN0] = self.calcSeqCB(v_ref, t_ref, Mi, setup)
        elif setup['Par']['PWM']['type'] == "OPP":
            [xs, xsh, s, c, x, xN0] = self.calcSeqOPP(v_ref, t_ref, Mi, setup)
        else:
            [xs, xsh, s, c, x, xN0] = self.calcSeqCB(v_ref, t_ref, Mi, setup)

        # ==============================================================================
        # Return
        # ==============================================================================
        return [xs, xsh, s, c, x, xN0]

    ###################################################################################################################
    # Temporal Output
    ###################################################################################################################
    def calcTime(self, s, e, t, Mi, mdl, t0, t1, setup):
        # ==============================================================================
        # Description
        # ==============================================================================
        """
        This function calculates the temporal outputs of the B4 bridge (AC and DC side).

        Input:
        1) s:       switching function
        2) e:       induced voltage
        3) t:       reference time (sec)
        4) Mi:      modulation index (0 ... 4/pi)
        5) mdl:     transfer functions
        6) t0:      start time (sample)
        7) t1:      end time (sample)
        8) setup:   all setup variables

        Output:
        1) outAc:   outputs ac side
        2) outDc:   outputs dc side
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
        # Inverter Output
        v_a0 = 0.5 * s['A'] * self.Vdc
        v_b0 = 0.5 * s['B'] * self.Vdc
        v_ab = v_a0 - v_b0

        # Filter Output
        if setup['Top']['outFilter'] == 0:
            v_out = v_ab
        else:
            _, v_out, _, = signal.lsim(mdl['SS']['Out'], v_ab, t)

        # Load
        v_L = v_out - Mi * e['A']

        # Current
        if setup['Top']['wave'] == "con":
            _, i_a, _, = signal.lsim(mdl['SS']['Load'], v_L, t)
            i_a = i_a[t0:t1]
        else:
            _, i_a, _, = signal.lsim(mdl['SS']['Load'], (v_L - np.mean(v_L)), t)
            i_a = i_a[t0:t1]
            i_a = i_a - np.mean(i_a)

        # ------------------------------------------
        # DC Side
        # ------------------------------------------
        # Inverter Input
        i_dc = i_a / 2 * (s['A'][t0:t1] - s['B'][t0:t1])

        # DC-Link
        i_c = np.mean(i_dc) - i_dc
        _, v_dc, _, = signal.lsim(mdl['SS']['DC'], i_c, t[t0:t1])

        # Filter Input
        if setup['Top']['inpFilter'] == 0:
            v_in = v_dc
        else:
            _, v_in, _, = signal.lsim(mdl['SS']['Inp'], (v_dc - self.Vdc), t[t0:t1])
            v_in = v_in + self.Vdc

        # ==============================================================================
        # Post-Processing
        # ==============================================================================
        # ------------------------------------------
        # AC Side
        # ------------------------------------------
        outAc['v_a0'] = v_a0[t0:t1]
        outAc['v_b0'] = v_a0[t0:t1]
        outAc['v_L'] = v_L[t0:t1]
        outAc['v_a_out'] = v_out[t0:t1]
        outAc['v_a'] = v_ab[t0:t1]
        outAc['i_a'] = i_a

        # ------------------------------------------
        # DC Side
        # ------------------------------------------
        outDc['v_in'] = v_in
        outDc['v_dc'] = v_dc
        outDc['i_dc'] = i_dc
        outDc['i_c'] = i_c

        # ==============================================================================
        # Return
        # ==============================================================================
        return [outAc, outDc]

    ###################################################################################################################
    # Analytical distortion
    ###################################################################################################################
    def calcDist(self, _1, _2, Mi, L, _3, _4):
        # ==============================================================================
        # Description
        # ==============================================================================
        """
        This function calculates the analytical distortion of the B4 bridge.

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
        # AC Side (Line-to-Neutral)
        # ------------------------------------------
        V_ab_eff = self.Vdc * np.sqrt(2 / np.pi * Mi) / 2
        V_ab_v1_eff = 1 / np.sqrt(2) * self.Vdc * Mi / 2
        V_ab_thd = self.Vdc * np.sqrt(1 - np.pi / 4 * Mi) / 2
        I_a_thd = 1 / np.sqrt(48) * self.Vdc * self.Ts / L * np.sqrt(3 / 8 * Mi ** 4 - 8 / (3 * np.pi) * Mi ** 3 + 0.5 * Mi ** 2)

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
        outAc['V_a_eff'] = V_ab_eff
        outAc['V_a_v1_eff'] = V_ab_v1_eff
        outAc['V_a_thd'] = V_ab_thd
        outAc['I_a_eff'] = np.nan
        outAc['I_a_v1_eff'] = np.nan
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
    def calcRef(self, E, phiE, phiV, setup):
        # ==============================================================================
        # Description
        # ==============================================================================
        """
        This function calculates the reference voltage and back EMF functions based on the
        B4 topology and the given parameters.

        Input:
        1) E:       amplitude of the back emf (V)
        2) phiE:    angle of the back emf (rad)
        3) v_a:     load angle of the output (rad)
        4) setup:   file including all setup parameters

        Output:
        1) v_ref:   reference voltage for given load scenario (V)
        2) e_ref:   reference back emf for given load scenario (V)
        """

        # ==============================================================================
        # Initialisation
        # ==============================================================================
        v_ref = {}
        e_ref = {}

        # ==============================================================================
        # Calculation
        # ==============================================================================
        # ------------------------------------------
        # Time
        # ------------------------------------------
        t = np.linspace(0, self.K / self.fel, self.K * self.N + 1)

        # ------------------------------------------
        # Reference
        # ------------------------------------------
        v_ref['A'] = +(self.Vdc / 2) * self.Mi * genWave(t, self.fel, phiV, setup)
        v_ref['B'] = -(self.Vdc / 2) * self.Mi * genWave(t, self.fel, phiV, setup)
        e_ref['A'] = E * genWave(t, self.fel, phiE, setup)
        e_ref['B'] = E * genWave(t, self.fel, phiE, setup)

        # ==============================================================================
        # Return
        # ==============================================================================
        return [v_ref, e_ref]

#######################################################################################################################
# References
#######################################################################################################################
