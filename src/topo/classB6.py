#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         classB6
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
This class initialises an object of the B6 (three-phase full-bridge) circuit.
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
from src.pwm.genSwSeq import genSwSeq
from src.pwm.svPWM import svPWM

# ==============================================================================
# External
# ==============================================================================
import numpy as np
import pandas as pd
from scipy import signal
from scipy.fft import fft


#######################################################################################################################
# Class
#######################################################################################################################
class classB6:
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
        self.Mi_max = 2 / np.sqrt(3)
        self.id1 = ['A', 'B', 'C']
        self.id2 = ['S1', 'S2', 'S3', 'S4', 'S5', 'S6']
        self.id3 = ['A', 'A', 'B', 'B', 'C', 'C']
        self.id4 = ['i_a', 'i_a', 'i_b', 'i_b', 'i_c', 'i_c']
        self.id5 = ['HS', 'LS', 'HS', 'LS', 'HS', 'LS']
        self.id6 = ['T1', 'T2', 'T3', 'T4', 'T5', 'T6']
        self.id7 = ['D1', 'D2', 'D3', 'D4', 'D5', 'D6']
        self.id8 = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6']
        self.id9 = [1, 1, 1, 1, 1, 1]
        self.name = 'B6'

    ###################################################################################################################
    # Init Data
    ###################################################################################################################
    def initData(self):
        # ==============================================================================
        # Description
        # ==============================================================================
        """
        This function initialises the data structure of the B6 bridge. This method is
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
        data['elec']['sw']['S5'] = pd.DataFrame(columns=['i_T', 'v_T', 'i_D', 'v_D'])
        data['elec']['sw']['S6'] = pd.DataFrame(columns=['i_T', 'v_T', 'i_D', 'v_D'])

        # ------------------------------------------
        # Losses
        # ------------------------------------------
        data['loss']['sw'] = {}
        data['loss']['sw']['S1'] = pd.DataFrame(columns=['p_T_c', 'p_T_s', 'p_T', 'p_D_c', 'p_D_s', 'p_D', 'p_L'])
        data['loss']['sw']['S2'] = pd.DataFrame(columns=['p_T_c', 'p_T_s', 'p_T', 'p_D_c', 'p_D_s', 'p_D', 'p_L'])
        data['loss']['sw']['S3'] = pd.DataFrame(columns=['p_T_c', 'p_T_s', 'p_T', 'p_D_c', 'p_D_s', 'p_D', 'p_L'])
        data['loss']['sw']['S4'] = pd.DataFrame(columns=['p_T_c', 'p_T_s', 'p_T', 'p_D_c', 'p_D_s', 'p_D', 'p_L'])
        data['loss']['sw']['S5'] = pd.DataFrame(columns=['p_T_c', 'p_T_s', 'p_T', 'p_D_c', 'p_D_s', 'p_D', 'p_L'])
        data['loss']['sw']['S6'] = pd.DataFrame(columns=['p_T_c', 'p_T_s', 'p_T', 'p_D_c', 'p_D_s', 'p_D', 'p_L'])

        # ------------------------------------------
        # Thermal
        # ------------------------------------------
        data['ther']['sw'] = pd.DataFrame(columns=['T1', 'T2', 'T3', 'T4', 'T5', 'T6',
                                                   'D1', 'D2', 'D3', 'D4', 'D5', 'D6',
                                                   'C1', 'C2', 'C3', 'C4', 'C5', 'C6'])

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
        This function initialises the output files of the B6 bridge

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
        timeSw = pd.DataFrame(columns=['t', 'v_a_ref', 'v_b_ref', 'v_c_ref', 'e_a', 'e_b', 'e_c', 'xAs', 'xBs', 'xCs',
                                       'xAsh', 'xBsh', 'xCsh', 'sA', 'sB', 'sC', 'xA', 'xB', 'xC', 'n0', 'c'])
        freqSw = pd.DataFrame(columns=['Sa', 'Xas'])
        freqAc = pd.DataFrame(columns=['I_a', 'V_a', 'V_a0'])
        freqDc = pd.DataFrame(columns=['I_dc', 'I_d_p', 'I_d_m', 'V_dc'])
        distAc['num'] = pd.DataFrame(data=np.ones((self.W, 6)), columns=['V_a_eff', 'V_a_v1_eff', 'V_a_thd', 'I_a_eff', 'I_a_v1_eff', 'I_a_thd'])
        distAc['ana'] = pd.DataFrame(data=np.ones((self.W, 6)), columns=['V_a_eff', 'V_a_v1_eff', 'V_a_thd', 'I_a_eff', 'I_a_v1_eff', 'I_a_thd'])
        distDc['num'] = pd.DataFrame(data=np.ones((self.W, 6)), columns=['V_dc_eff', 'V_dc_v1_eff', 'V_dc_thd', 'I_dc_eff', 'I_dc_v1_eff', 'I_dc_thd'])
        distDc['ana'] = pd.DataFrame(data=np.ones((self.W, 6)), columns=['V_dc_eff', 'V_dc_v1_eff', 'V_dc_thd', 'I_dc_eff', 'I_dc_v1_eff', 'I_dc_thd'])

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
        This function summarizes the output structure of the B6 bridge.

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
        timeSw['v_c_ref'] = v_ref['C'][t0:t1]
        timeSw['e_a'] = e_ref['A'][t0:t1]
        timeSw['e_b'] = e_ref['B'][t0:t1]
        timeSw['e_c'] = e_ref['C'][t0:t1]
        timeSw['sA'] = s['A'][t0:t1]
        timeSw['sB'] = s['B'][t0:t1]
        timeSw['sC'] = s['C'][t0:t1]
        timeSw['c'] = c[t0:t1]
        timeSw['xAs'] = xs['A'][t0:t1]
        timeSw['xBs'] = xs['B'][t0:t1]
        timeSw['xCs'] = xs['C'][t0:t1]
        timeSw['xAsh'] = xsh['A'][t0:t1]
        timeSw['xBsh'] = xsh['B'][t0:t1]
        timeSw['xCsh'] = xsh['C'][t0:t1]
        timeSw['xA'] = x['A'][t0:t1]
        timeSw['xB'] = x['B'][t0:t1]
        timeSw['xC'] = x['C'][t0:t1]
        timeSw['n0'] = xN0[t0:t1]

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
        s = {}
        x = {}

        # ==============================================================================
        # Calculation
        # ==============================================================================
        # ------------------------------------------
        # Reference
        # ------------------------------------------
        x0 = {'A': Mi * v_ref['A'] / np.max(np.abs(v_ref['A'])),
              'B': Mi * v_ref['B'] / np.max(np.abs(v_ref['B'])),
              'C': Mi * v_ref['C'] / np.max(np.abs(v_ref['C']))}

        # ------------------------------------------
        # Carrier
        # ------------------------------------------
        c = np.zeros(np.size(t_ref))

        # ------------------------------------------
        # Sequence
        # ------------------------------------------
        for i in range(0, len(self.id1)):
            s[self.id1[i]] = signal.square(2 * np.pi * self.fel * t_ref - i * (np.pi * 2) / 3, duty=Mi / 2)

        # ------------------------------------------
        # Neutral
        # ------------------------------------------
        xN0 = 1 / 3 * (s['A'] + s['B'] + s['C'])
        for i in range(0, len(self.id1)):
            x[self.id1[i]] = x0[self.id1[i]] + xN0

        # ------------------------------------------
        # Output
        # ------------------------------------------
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
        s = {}
        xs = {}
        xsh = {}

        # ------------------------------------------
        # Parameters
        # ------------------------------------------
        tmin = int(self.tmin / (t_ref[1] - t_ref[0]))

        # ==============================================================================
        # Pre-processing
        # ==============================================================================
        # ------------------------------------------
        # Reference
        # ------------------------------------------
        x0 = {'A': Mi * v_ref['A'] / np.max(np.abs(v_ref['A'])),
              'B': Mi * v_ref['B'] / np.max(np.abs(v_ref['B'])),
              'C': Mi * v_ref['C'] / np.max(np.abs(v_ref['C']))}
        xN0 = np.zeros(np.size(x0['A']))
        xAll = np.vstack((x0['A'], x0['B'], x0['C'])).transpose()

        # ------------------------------------------
        # Clark Transform
        # ------------------------------------------
        phi = np.arcsin(v_ref['A'][0] / np.max(v_ref['A']))

        # ------------------------------------------
        # Zero-Sequence
        # ------------------------------------------
        # SPWM
        if setup['Par']['PWM']['zero'] == "SPWM":
            xN0 = np.zeros(np.size(x0['A']))

        # SVPWM
        if setup['Par']['PWM']['zero'] == "SVPWM":
            xN0 = 1 / 4 * Mi * signal.sawtooth(3 * 2 * np.pi * self.fel * (t_ref - (0.25 - phi / (2 * np.pi)) / self.fel), 0.5)

        # THIPWM1/4
        if setup['Par']['PWM']['zero'] == "THIPWM4":
            xN0 = 1 / 4 * Mi * np.sin(3 * 2 * np.pi * self.fel * (t_ref + (phi / (2 * np.pi)) / self.fel))

        # THIPWM1/6
        if setup['Par']['PWM']['zero'] == "THIPWM6":
            xN0 = 1 / 4 * Mi * np.sin(3 * 2 * np.pi * self.fel * (t_ref + (phi / (2 * np.pi)) / self.fel))

        # DPWM0
        if setup['Par']['PWM']['zero'] == "DPWM0":
            xAll_s = np.roll(xAll, shift=-int(len(xN0) / 12), axis=0)
            id2 = np.argsort(abs(xAll_s), axis=1)
            for i in range(0, len(xN0)):
                xN0[i] = np.sign(xAll[i, id2[i, 2]]) - xAll[i, id2[i, 2]]

        # DPWM1
        if setup['Par']['PWM']['zero'] == "DPWM1":
            xAll_s = np.roll(xAll, shift=0, axis=0)
            id2 = np.argsort(abs(xAll_s), axis=1)
            for i in range(0, len(xN0)):
                xN0[i] = np.sign(xAll[i, id2[i, 2]]) - xAll[i, id2[i, 2]]

        # DPWM2
        if setup['Par']['PWM']['zero'] == "DPWM2":
            xAll_s = np.roll(xAll, shift=int(len(xN0) / 12), axis=0)
            id2 = np.argsort(abs(xAll_s), axis=1)
            for i in range(0, len(xN0)):
                xN0[i] = np.sign(xAll[i, id2[i, 2]]) - xAll[i, id2[i, 2]]

        # DPWM3
        if setup['Par']['PWM']['zero'] == "DPWM3":
            id2 = np.argsort(abs(xAll), axis=1)
            for i in range(0, len(xN0)):
                xN0[i] = np.sign(xAll[i, id2[i, 1]]) - xAll[i, id2[i, 1]]

        # DPWMMIN
        if setup['Par']['PWM']['zero'] == "DPWMMIN":
            xN0 = -1 - np.min(xAll, axis=1)

        # DPWMMAX
        if setup['Par']['PWM']['zero'] == "DPWMMAX":
            xN0 = 1 - np.max(xAll, axis=1)

        # ------------------------------------------
        # Line-to-Line References
        # ------------------------------------------
        x['A'] = x0['A'] + xN0
        x['B'] = x0['B'] + xN0
        x['C'] = x0['C'] + xN0

        # ==============================================================================
        # Calculation
        # ==============================================================================
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
            c = signal.sawtooth(2 * np.pi * self.fs * (t_ref - 0.5 / self.fs), 0.5)
        c = (2 * (c - min(c)) / (max(c) - min(c))) - 1

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
        for i in range(0, len(self.id1)):
            s[self.id1[i]] = cbInter(xs[self.id1[i]], c, Mi, tmin)

        # ------------------------------------------
        # Dead-Time
        # ------------------------------------------
        if setup['Par']['PWM']['td'] > 0:
            for i in range(0, len(self.id1)):
                s[self.id1[i]] = self.calcDead(s[self.id1[i]], t_ref, Mi)

        # ==============================================================================
        # Return
        # ==============================================================================
        return [xs, xsh, s, c, x, xN0]

    ###################################################################################################################
    # Optimal Pulse Patterns
    ###################################################################################################################
    def calcSeqSV(self, v_ref, t_ref, Mi, setup):
        # ==============================================================================
        # Description
        # ==============================================================================
        """
        This function calculates the switching time assuming space vector modulation.

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
        # Parameters
        # ------------------------------------------
        N1 = int(len(t_ref) / (self.q * self.K))
        if setup['Par']['PWM']['upd'] == "SE":
            Ns = self.q
            Terr = -N1 / 2
        else:
            Ns = 2 * self.q
            Terr = -N1 / 4

        # ------------------------------------------
        # Variables
        # ------------------------------------------
        # Empty
        x = {}
        s = {}
        xs = {}
        xsh = {}
        x0 = {}
        mS = {}

        # Zeros
        s['A'] = np.zeros(np.size(t_ref))
        s['B'] = np.zeros(np.size(t_ref))
        s['C'] = np.zeros(np.size(t_ref))
        xs['A'] = np.zeros(np.size(t_ref))
        xs['B'] = np.zeros(np.size(t_ref))
        xs['C'] = np.zeros(np.size(t_ref))
        ts = np.linspace(0, 2, N1)
        ss = np.zeros(np.size(t_ref))
        c = np.zeros(np.size(t_ref))
        xN0 = np.zeros(np.size(t_ref))
        rr = np.zeros((Ns * self.K, 1))
        t1 = np.zeros((Ns * self.K, 1))
        t2 = np.zeros((Ns * self.K, 1))
        t0 = np.zeros((Ns * self.K, 1))
        t7 = np.zeros((Ns * self.K, 1))

        # Mapping
        mS['A'] = [-1, +1, +1, -1, -1, -1, +1, +1]
        mS['B'] = [-1, -1, +1, +1, +1, -1, -1, +1]
        mS['C'] = [-1, -1, -1, -1, +1, +1, +1, +1]

        # ==============================================================================
        # Calculation
        # ==============================================================================
        # ------------------------------------------
        # Reference
        # ------------------------------------------
        x0['A'] = Mi * v_ref['A'] / np.max(v_ref['A'])
        x0['B'] = Mi * v_ref['B'] / np.max(v_ref['B'])
        x0['C'] = Mi * v_ref['C'] / np.max(v_ref['C'])

        # ------------------------------------------
        # Clark Transform
        # ------------------------------------------
        v_ref['alpha'] = 2 / 3 * (v_ref['A'] - 0.5 * v_ref['B'] - 0.5 * v_ref['C'])
        v_ref['beta'] = 2 / 3 * (np.sqrt(3) / 2 * v_ref['B'] - np.sqrt(3) / 2 * v_ref['C'])
        phi = np.arctan2(v_ref['beta'], v_ref['alpha'])

        # ------------------------------------------
        # Define Switching Sequence
        # ------------------------------------------
        [seq, k] = genSwSeq(setup)

        # ------------------------------------------
        # Zero-Sequence
        # ------------------------------------------
        for i in range(0, len(t_ref)):
            alpha = i * self.K / len(t_ref) * 2 * np.pi + phi[0] + 2 * np.pi
            [d0, d1, d2, d7, _] = svPWM(k, alpha, Mi)
            xN0[i] = (-d0 - d1 / 3 + d2 / 3 + d7)

        # ------------------------------------------
        # Line-to-Line References
        # ------------------------------------------
        x['A'] = x0['A'] + xN0
        x['B'] = x0['B'] + xN0
        x['C'] = x0['C'] + xN0

        # ==============================================================================
        # Calculation
        # ==============================================================================
        # ------------------------------------------
        # Determine Sector
        # ------------------------------------------
        for i in range(0, Ns * self.K):
            alpha = i / Ns * 2 * np.pi + phi[0] + 2 * np.pi
            alpha = alpha + Terr / len(t_ref) * self.K * (2 * np.pi)
            [t0[i], t1[i], t2[i], t7[i], rr[i]] = svPWM(k, alpha, Mi)

        # ------------------------------------------
        # Switching times
        # ------------------------------------------
        if setup['Par']['PWM']['upd'] == "SE":
            st = np.hstack((t0, t0 + t1, 1 - t7, np.ones((Ns * self.K, 1)), 1 + t7, 1 + t7 + t2, 2 - t0, 2 * np.ones((Ns * self.K, 1))))
        else:
            st1 = np.hstack((t0, t0 + t1, 1 - t7, np.ones((Ns * self.K, 1))))
            st2 = np.roll(np.hstack((1 + t7, 1 + t7 + t2, 2 - t0, 2 * np.ones((Ns * self.K, 1)))), -1, axis=0)
            st = np.hstack((st1, st2))
            st = st[::2]
            rr = rr[::2]

        # ------------------------------------------
        # Switching states
        # ------------------------------------------
        for i in range(0, self.q * self.K):
            j = 0
            for ii in range(0, N1):
                if st[i, j] > ts[ii]:
                    ss[ii + N1 * i] = seq[int(rr[i] - 1)][j]
                else:
                    j = j + 1
                    ss[ii + N1 * i] = ss[ii + N1 * i - 1]

        # ------------------------------------------
        # Sampling
        # ------------------------------------------
        for i in range(0, len(self.id1)):
            for ii in range(0, len(x[self.id1[i]])):
                if (ii % int(len(t_ref) / (Ns * self.K))) == 0:
                    xs[self.id1[i]][ii] = x[self.id1[i]][ii]
                else:
                    xs[self.id1[i]][ii] = xs[self.id1[i]][ii - 1]

        # ------------------------------------------
        # Shifting
        # ------------------------------------------
        for i in range(0, len(self.id1)):
            if setup['Par']['PWM']['upd'] == "SE":
                xsh[self.id1[i]] = np.roll(x0[self.id1[i]], int(len(xs[self.id1[i]]) * self.fel / self.fs))
            else:
                xsh[self.id1[i]] = np.roll(x0[self.id1[i]], int(len(xs[self.id1[i]]) * self.fel / self.fs / 2))

        # ------------------------------------------
        # Mapping
        # ------------------------------------------
        for i in range(0, len(self.id1)):
            for ii in range(0, len(ss)):
                if Mi != 0:
                    s[self.id1[i]][ii] = mS[self.id1[i]][int(ss[ii])]

        # ------------------------------------------
        # Dead time
        # ------------------------------------------
        for i in range(0, len(self.id1)):
            s[self.id1[i]] = self.calcDead(s[self.id1[i]], t_ref, Mi)

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
        s = {}
        xs = {}
        xsh = {}
        x0 = {}
        ang_total = []
        val_total = []

        # Zeros
        s['A'] = np.zeros(np.size(t_ref))
        s['B'] = np.zeros(np.size(t_ref))
        s['C'] = np.zeros(np.size(t_ref))
        xs['A'] = np.zeros(np.size(t_ref))
        xs['B'] = np.zeros(np.size(t_ref))
        xs['C'] = np.zeros(np.size(t_ref))
        c = np.zeros(np.size(t_ref))
        ss = (-1) * np.ones(np.size(t_ref))
        xN0 = np.zeros(np.size(t_ref))

        # ------------------------------------------
        # Parameters
        # ------------------------------------------
        kmax = 10 * self.q
        ang = np.linspace(0, np.pi * 2 * self.N, np.size(t_ref))
        if setup['Par']['PWM']['upd'] == "SE":
            Ns = self.q
        else:
            Ns = 2 * self.q

        # ==============================================================================
        # Calculation
        # ==============================================================================
        # ------------------------------------------
        # Reference
        # ------------------------------------------
        x0['A'] = Mi * v_ref['A'] / np.max(v_ref['A'])
        x0['B'] = Mi * v_ref['B'] / np.max(v_ref['B'])
        x0['C'] = Mi * v_ref['C'] / np.max(v_ref['C'])
        x['A'] = x0['A'] + xN0
        x['B'] = x0['B'] + xN0
        x['C'] = x0['C'] + xN0

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
        for i in range(1, len(ss)):
            if c[i] == 0:
                ss[i] = ss[i - 1]
            else:
                ss[i] = c[i]

        # Direction
        for i in range(1, len(ss)):
            if v_ref['A'][i] >= 0:
                ss[i] = ss[i] * (+1)
            else:
                ss[i] = ss[i] * (-1)
        ss = -ss

        # ------------------------------------------
        # Three phase
        # ------------------------------------------
        for i in range(0, len(self.id1)):
            s[self.id1[i]] = np.roll(ss, int(np.floor(i * 120 / 360 / self.N * len(ss))))

        # ------------------------------------------
        # Dead time
        # ------------------------------------------
        for i in range(0, len(self.id1)):
            s[self.id1[i]] = self.calcDead(s[self.id1[i]], t_ref, Mi)

        # ------------------------------------------
        # Sampling
        # ------------------------------------------
        for i in range(0, len(self.id1)):
            for ii in range(0, len(x[self.id1[i]])):
                if (ii % int(len(t_ref) / (Ns * self.N))) == 0:
                    xs[self.id1[i]][ii] = x[self.id1[i]][ii]
                else:
                    xs[self.id1[i]][ii] = xs[self.id1[i]][ii - 1]

        # ------------------------------------------
        # Shifting
        # ------------------------------------------
        for i in range(0, len(self.id1)):
            if setup['Par']['PWM']['upd'] == "SE":
                xsh[self.id1[i]] = np.roll(x0[self.id1[i]], int(len(xs[self.id1[i]]) * self.fel / self.fs))
            else:
                xsh[self.id1[i]] = np.roll(x0[self.id1[i]], int(len(xs[self.id1[i]]) * self.fel / self.fs / 2))

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
        PWM methods are carrier-based, fundamental frequency, space vector modulation,
        and optimal pulse patterns.

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
        elif setup['Par']['PWM']['type'] == "SV":
            [xs, xsh, s, c, x, xN0] = self.calcSeqSV(v_ref, t_ref, Mi, setup)
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
    def calcTime(self, s, e, t, Mi, mdl, t0, t1, init, avg, setup):
        # ==============================================================================
        # Description
        # ==============================================================================
        """
        This function calculates the temporal outputs of the B6 bridge (AC and DC side).

        Input:
        1) s:       switching function
        2) e:       induced voltage
        3) t:       reference time (sec)
        4) Mi:      modulation index (0 ... 4/pi)
        5) mdl:     transfer functions
        6) t0:      start time (sample)
        7) t1:      end time (sample)
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
        # Variables
        # ------------------------------------------
        v0 = {}
        v = {}
        v_out = {}
        v_L = {}
        i = {}

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
            init = {'inp': 0, 'out': [0, 0, 0], 'dc': 0, 'load': 0}

        # ------------------------------------------
        # Parameters
        # ------------------------------------------
        K1 = int((t1 - t0 - 1) * self.fel * (t[1] - t[0]))

        # ==============================================================================
        # Calculation
        # ==============================================================================
        # ------------------------------------------
        # AC Side
        # ------------------------------------------
        # Inverter Output
        for j in range(0, len(self.id1)):
            v0[self.id1[j]] = 0.5 * (s[self.id1[j]] - np.mean(s[self.id1[j]])) * self.Vdc
        v_n0 = 1 / 3 * (v0['A'] + v0['B'] + v0['C'])

        # Phase Voltages
        for j in range(0, len(self.id1)):
            v[self.id1[j]] = v0[self.id1[j]] - v_n0

        # Filter Output
        for j in range(0, len(self.id1)):
            if setup['Top']['outFilter'] == 0:
                v_out[self.id1[j]] = v[self.id1[j]]
            else:
                _, v_out[self.id1[j]], _, = signal.lsim(mdl['SS']['Out'], v[self.id1[j]], t, X0=init['out'][j])

        # Load
        for j in range(0, len(self.id1)):
            v_L[self.id1[j]] = v_out[self.id1[j]] - Mi * e[self.id1[j]]

        # LL Current
        _, i_a, _, = signal.lsim(mdl['SS']['Load'], v_L['A'], t, X0=init['load'])
        i['A'] = i_a[t0:t1]
        i['B'] = np.roll(i_a[t0:t1], int(np.floor(120 / 360 / K1 * len(s['A'][t0:t1]))))
        i['C'] = np.roll(i_a[t0:t1], int(np.floor(240 / 360 / K1 * len(s['A'][t0:t1]))))

        # LN Current
        if setup['Top']['wave'] != 'con':
            for j in range(0, len(self.id1)):
                if avg == 1:
                    i[self.id1[j]] = i[self.id1[j]] - np.mean(i[self.id1[j]])

        # ------------------------------------------
        # DC Side
        # ------------------------------------------
        # Inverter Input
        i_dc = 1 / 2 * (s['A'][t0:t1] * i['A'] + s['B'][t0:t1] * i['B'] + s['C'][t0:t1] * i['C'])

        # DC-Link
        i_cap = np.mean(i_dc) - i_dc
        _, v_dc, _, = signal.lsim(mdl['SS']['DC'], i_cap, t[t0:t1], X0=init['dc'])
        v_dc = v_dc - np.mean(v_dc) + self.Vdc

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
        outAc['v_a0'] = v0['A'][t0:t1]
        outAc['v_a'] = v['A'][t0:t1]
        outAc['v_L_a'] = v_L['A'][t0:t1]
        outAc['v_a_out'] = v_out['A'][t0:t1]
        outAc['v_n0'] = v_n0[t0:t1]
        outAc['i_a'] = i['A']
        outAc['i_b'] = i['B']
        outAc['i_c'] = i['C']

        # ------------------------------------------
        # DC Side
        # ------------------------------------------
        outDc['v_in'] = v_in
        outDc['v_dc'] = v_dc
        outDc['i_dc'] = i_dc
        outDc['i_c'] = i_cap

        # ------------------------------------------
        # Init Conditions
        # ------------------------------------------
        outIn['inp'] = v_in[-1]
        outIn['out'] = [v_out['A'][-1], v_out['B'][-1], v_out['C'][-1]]
        outIn['dc'] = v_dc[-1]
        outIn['load'] = i_a[-1]

        # ==============================================================================
        # Return
        # ==============================================================================
        return [outAc, outDc, outIn]

    ###################################################################################################################
    # Analytical distortion
    ###################################################################################################################
    def calcDist(self, i_a, v_a, Mi, L, Z, setup):
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
        I_a_thd = []
        HDF = np.zeros(np.size(Mi))

        # ==============================================================================
        # Pre-processing
        # ==============================================================================
        Y = fft(v_a)
        phiV1 = np.angle(Y)[self.K]
        Y = fft(i_a)
        phiI1 = np.angle(Y)[self.K]
        phi = phiV1 - phiI1

        # ==============================================================================
        # Calculation
        # ==============================================================================
        # ------------------------------------------
        # AC Side
        # ------------------------------------------
        # Voltages
        V_a_eff = self.Vdc * np.sqrt(Mi / (np.sqrt(3) * np.pi))
        V_a_v1_eff = 1 / np.sqrt(2) * self.Vdc / 2 * Mi
        V_a_thd = self.Vdc / 2 * np.sqrt(1 - Mi * (np.sqrt(3) * np.pi) / 8)

        # Currents
        I_a_eff = V_a_eff / Z
        I_a_v1_eff = V_a_v1_eff / Z

        # Distortion
        if setup['Par']['PWM']['type'] == "FF":
            I_a_thd = 0.0417 * (self.Vdc / 2) * (1 / (2 * np.pi * self.fel * L)) * Mi
        else:
            # 0127 (SPWM)
            if setup['Par']['PWM']['seq'] == "0127" and setup['Par']['PWM']['zero'] == "SPWM":
                HDF = 3 / 2 * Mi ** 2 - 4 * np.sqrt(3) / np.pi * Mi ** 3 + 9 / 8 * Mi ** 4

            # 0127 (SVPWM)
            if setup['Par']['PWM']['seq'] == "0127" and setup['Par']['PWM']['zero'] == "SVPWM":
                HDF = 3 / 2 * Mi ** 2 - 4 * np.sqrt(3) / np.pi * Mi ** 3 + (
                        27 / 16 - 81 * np.sqrt(3) / (64 * np.pi)) * Mi ** 4

            # 0127 (THIPWM4)
            if setup['Par']['PWM']['seq'] == "0127" and setup['Par']['PWM']['zero'] == "THIPWM4":
                HDF = 3 / 2 * Mi ** 2 - 4 * np.sqrt(3) / np.pi * Mi ** 3 + 63 / 64 * Mi ** 4

            # 0127 (THIPWM6)
            if setup['Par']['PWM']['seq'] == "0127" and setup['Par']['PWM']['zero'] == "THIPWM4":
                HDF = 3 / 2 * Mi ** 2 - 4 * np.sqrt(3) / np.pi * Mi ** 3 + Mi ** 4

            # 0127 (DPWM0)
            if setup['Par']['PWM']['seq'] == "0127" and setup['Par']['PWM']['zero'] == "DPWM0":
                HDF_max = 6 * Mi ** 2 - (8 * np.sqrt(3) + 45) / (2 * np.pi) * Mi ** 3 + (
                        27 / 8 + 27 * np.sqrt(3) / (32 * np.pi)) * Mi ** 4
                HDF_min = 6 * Mi ** 2 + (45 - 62 * np.sqrt(3)) / (2 * np.pi) * Mi ** 3 + (
                        27 / 8 + 27 * np.sqrt(3) / (16 * np.pi)) * Mi ** 4
                HDF = 0.5 * (HDF_max + HDF_min)

            # 0127 (DPWM1)
            if setup['Par']['PWM']['seq'] == "0127" and setup['Par']['PWM']['zero'] == "DPWM1":
                HDF = 6 * Mi ** 2 - (8 * np.sqrt(3) + 45) / (2 * np.pi) * Mi ** 3 + (
                        27 / 8 + 27 * np.sqrt(3) / (32 * np.pi)) * Mi ** 4

            # 0127 (DPWM2)
            if setup['Par']['PWM']['seq'] == "0127" and setup['Par']['PWM']['zero'] == "DPWM2":
                HDF_max = 6 * Mi ** 2 - (8 * np.sqrt(3) + 45) / (2 * np.pi) * Mi ** 3 + (
                        27 / 8 + 27 * np.sqrt(3) / (32 * np.pi)) * Mi ** 4
                HDF_min = 6 * Mi ** 2 + (45 - 62 * np.sqrt(3)) / (2 * np.pi) * Mi ** 3 + (
                        27 / 8 + 27 * np.sqrt(3) / (16 * np.pi)) * Mi ** 4
                HDF = 0.5 * (HDF_max + HDF_min)

            # 0127 (DPWM3)
            if setup['Par']['PWM']['seq'] == "0127" and setup['Par']['PWM']['zero'] == "DPWM3":
                HDF = 6 * Mi ** 2 + (45 - 62 * np.sqrt(3)) / (2 * np.pi) * Mi ** 3 + (
                        27 / 8 + 27 * np.sqrt(3) / (16 * np.pi)) * Mi ** 4

            # 0127 (DPWMMAX)
            if setup['Par']['PWM']['seq'] == "0127" and setup['Par']['PWM']['zero'] == "DPWMMAX":
                HDF_max = 6 * Mi ** 2 - (8 * np.sqrt(3) + 45) / (2 * np.pi) * Mi ** 3 + (
                        27 / 8 + 27 * np.sqrt(3) / (32 * np.pi)) * Mi ** 4
                HDF_min = 6 * Mi ** 2 + (45 - 62 * np.sqrt(3)) / (2 * np.pi) * Mi ** 3 + (
                        27 / 8 + 27 * np.sqrt(3) / (16 * np.pi)) * Mi ** 4
                HDF = 0.5 * (HDF_max + HDF_min)

            # 0127 (DPWMMIN)
            if setup['Par']['PWM']['seq'] == "0127" and setup['Par']['PWM']['zero'] == "DPWMMIN":
                HDF_max = 6 * Mi ** 2 - (8 * np.sqrt(3) + 45) / (2 * np.pi) * Mi ** 3 + (
                        27 / 8 + 27 * np.sqrt(3) / (32 * np.pi)) * Mi ** 4
                HDF_min = 6 * Mi ** 2 + (45 - 62 * np.sqrt(3)) / (2 * np.pi) * Mi ** 3 + (
                        27 / 8 + 27 * np.sqrt(3) / (16 * np.pi)) * Mi ** 4
                HDF = 0.5 * (HDF_max + HDF_min)

        # ------------------------------------------
        # DC Side
        # ------------------------------------------
        # Voltages
        V_dc_thd = np.sqrt(Mi ** 2 * (np.sqrt(3) / (4 * np.pi) + np.cos(phi) ** 2 * (np.sqrt(3) / np.pi - 6 / 16 * Mi))) * np.sqrt(I_a_v1_eff)

        # Currents
        I_dc_eff = I_a_v1_eff * np.sqrt(2 * np.sqrt(3) / np.pi * Mi * (1 / 4 + np.cos(phi) ** 2))
        I_dc_v1_eff = 3 / 4 * np.sqrt(2) * I_a_v1_eff * np.cos(phi)
        I_dc_thd = np.sqrt(2 * Mi * (np.sqrt(3) / (4 * np.pi) + np.cos(phi) ** 2 * (np.sqrt(3) / np.pi - 9 / 16 * Mi))) * I_a_v1_eff

        # ==============================================================================
        # Post-Processing
        # ==============================================================================
        # ------------------------------------------
        # Denormalization
        # ------------------------------------------
        if setup['Par']['PWM']['type'] != "FF":
            I_a_thd = self.Vdc / (24 * L * self.fs) * np.sqrt(HDF)

        # ------------------------------------------
        # AC Side
        # ------------------------------------------
        outAc['V_a_eff'] = V_a_eff
        outAc['V_a_v1_eff'] = V_a_v1_eff
        outAc['V_a_thd'] = V_a_thd
        outAc['I_a_eff'] = I_a_eff
        outAc['I_a_v1_eff'] = I_a_v1_eff
        outAc['I_a_thd'] = I_a_thd

        # ------------------------------------------
        # DC Side
        # ------------------------------------------
        outDc['V_dc_eff'] = np.nan
        outDc['V_dc_v1_eff'] = np.nan
        outDc['V_dc_thd'] = np.nan
        outDc['I_dc_eff'] = I_dc_eff
        outDc['I_dc_v1_eff'] = I_dc_v1_eff
        outDc['I_dc_thd'] = I_dc_thd

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
        3) i_ref:   reference current for given load scenario (A)
        """

        # ==============================================================================
        # Initialisation
        # ==============================================================================
        v_ref = {}
        e_ref = {}
        i_ref = {}

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
        for i in range(0, len(self.id1)):
            v_ref[self.id1[i]] = (self.Vdc / 2) * self.Mi * genWave(t, self.fel, phiV - i * 2 / 3 * np.pi, setup)
            e_ref[self.id1[i]] = E * genWave(t, self.fel, phiE - i * 2 / 3 * np.pi, setup)
            i_ref[self.id1[i]] = setup['Dat']['stat']['Io'] * np.sqrt(2) * genWave(t, self.fel, setup['Dat']['stat']['PhiVI'] - i * 2 / 3 * np.pi, setup)

        # ==============================================================================
        # Return
        # ==============================================================================
        return [v_ref, e_ref, i_ref]

#######################################################################################################################
# References
#######################################################################################################################
