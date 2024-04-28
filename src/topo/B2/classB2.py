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
Inputs:     1)
            2)
            N)
"""

#######################################################################################################################
# Import libs
#######################################################################################################################
# ==============================================================================
# Internal
# ==============================================================================
from src.general.helpFnc import cbInter, con2dis, deadTime
from src.pwm.oppPWM import oppPWM

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
class classB2:
    ###################################################################################################################
    # Constructor
    ###################################################################################################################
    def __init__(self, fel, fs, fsim, td, tmin, cyc, W, Mi, Vdc, Tc, Tj):
        self.fel = fel
        self.fs = fs
        self.fsim = fsim
        self.td = td
        self.tmin = tmin
        self.cyc = cyc
        self.W = W
        self.Mi = Mi
        self.Vdc = Vdc
        self.Tc = Tc
        self.Tj = Tj
        self.Nsim = int(np.ceil(fsim/fel))
        self.Npwm = int(np.ceil(fs/fel))
        self.q = int(fs / fel)
        self.N = int(fel / fsim)
        self.Ts = 1 / fs
        self.Tel = 1 / fel

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
        1) timeElec:    time domain electrical
        1) timeLoss:    time domain losses
        1) timeTher:    time domain thermal
        1) freqSw:      frequency domain switching function
        1) freqDc:      frequency domain dc
        1) freqAc:      frequency domain ac
        1) distAc:      distortion ac
        1) distDc:      distortion dc
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
        timeSw = pd.DataFrame(columns=['t', 'v_ref', 'e', 'xs', 'xsh', 's', 'c'])
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
        """

        # ==============================================================================
        # Initialisation
        # ==============================================================================
        c = np.zeros(np.size(t_ref))

        # ==============================================================================
        # Calculation
        # ==============================================================================
        x = Mi * v_ref / np.max(v_ref)
        s = signal.square(2 * np.pi * self.fel * t_ref, duty=0.5)
        xs = x
        xsh = x

        # ==============================================================================
        # Return
        # ==============================================================================
        return [xs, xsh, s, c]

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
        return [xs, xsh, s, c]

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

        # ------------------------------------------
        # Optimal Angles
        # ------------------------------------------
        # Fundamental angle (0 ... 2pi)
        [ang_fun, val_fun, _] = oppPWM(kmax, self.q, Mi / 4 * np.pi, 4, setup['Top'])

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
        return [xs, xsh, s, c]

    ###################################################################################################################
    # Temporal Output
    ###################################################################################################################
    def calcTime(self, s, e, t, Mi, mdl, t0, t1, setup):
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
        v_a0 = 0.5 * s * self.Vdc

        # Filter Output
        if setup['Top']['outFilter'] == 0:
            v_L = v_a0
        else:
            _, v_L, _, = signal.lsim(mdl['SS']['Out'], v_a0, t)

        # Load
        v_a = v_L - Mi * e

        # Current
        if setup['Top']['wave'] == "con":
            _, i_a, _, = signal.lsim(mdl['SS']['Load'], v_a, t)
            i_a = i_a[t0:t1]
        else:
            _, i_a, _, = signal.lsim(mdl['SS']['Load'], (v_a - np.mean(v_a)), t)
            i_a = i_a[t0:t1]
            i_a = i_a - np.mean(i_a)

        # ------------------------------------------
        # DC Side
        # ------------------------------------------
        # Inverter Input
        i_d_p = i_a * (1 + s[t0:t1]) / 2
        i_d_m = i_a * (1 - s[t0:t1]) / 2
        i_dc = i_d_p

        # DC-Link
        i_c = np.mean(i_d_p) - i_d_p
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
        outAc['v_L'] = v_L[t0:t1]
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

        # ==============================================================================
        # Return
        # ==============================================================================
        return [outAc, outDc]

    ###################################################################################################################
    # Analytical distortion
    ###################################################################################################################
    def calcDist(self, Mi, L, Z):
        # ==============================================================================
        # Description
        # ==============================================================================
        """
        This function calculates the analytical distortion of the B2 bridge.

        Input:
        1) Mi:      modulation index (0 ... 4/pi)
        2) L:       load inductance
        3) Z:       load impedance

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
        I_dc_thd = 1 / np.sqrt(48) * self.Vdc / 2 * self.Ts / L * np.sqrt(3 / 8 * Mi ** 4 - Mi ** 2 + 1)

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
        outDc['V_dc_v1_eff'] = 0
        outDc['V_dc_thd'] = 0
        outDc['I_dc_eff'] = 0
        outDc['I_dc_v1_eff'] = 0
        outDc['I_dc_thd'] = I_dc_thd

        # ==============================================================================
        # Return
        # ==============================================================================
        return [outAc, outDc]

    ###################################################################################################################
    # Analytical distortion
    ###################################################################################################################
    def calcDistNum(self, t, i_a, v_a, i_dc, v_dc):
        # ==============================================================================
        # Description
        # ==============================================================================
        """
        This function calculates the numerical distortion of the B2 bridge current and
        voltage waveforms.

        Input:
        1) t:       input time vector (sec)
        2) i_a:     load current (A)
        3) v_a:     load voltage (V)
        4) i_dc:    dc current (A)
        5) v_dc:    dc voltage (V)

        Output:
        1) outAc:   outputs distortion ac side
        2) outDc:   outputs distortion dc side
        """

        # ==============================================================================
        # Initialisation
        # ==============================================================================
        dt = t[1] - t[0]
        K = int(np.round((t[-1] - t[0]) / self.Tel))
        N = int(len(v_a))
        outAc = {}
        outDc = {}

        # ==============================================================================
        # Calculation
        # ==============================================================================
        # ------------------------------------------
        # AC Side
        # ------------------------------------------
        V_a_eff = np.sqrt(1 / self.Tel / K * np.sum(v_a ** 2 * dt))
        V_a_v1_eff = (1 / np.sqrt(2)) * 2 * np.abs(fft(v_a) / N)[K]
        V_a_thd = np.sqrt(V_a_eff ** 2 - V_a_v1_eff ** 2) / V_a_eff * self.Vdc / 2
        I_a_eff = np.sqrt(1 / self.Tel / K * np.sum(i_a ** 2 * dt))
        I_a_v1_eff = (1 / np.sqrt(2)) * 2 * np.abs(fft(i_a) / N)[K]
        I_a_thd = np.sqrt(I_a_eff ** 2 - I_a_v1_eff ** 2)

        # ------------------------------------------
        # DC Side
        # ------------------------------------------
        V_dc_eff = np.sqrt(1 / self.Tel / K * np.sum(v_dc ** 2 * dt))
        V_dc_v1_eff = np.abs(fft(v_dc) / N)[0]
        V_dc_thd = np.sqrt((np.sqrt(1 / self.Tel / K * np.sum((v_dc - self.Vdc) ** 2 * dt))) ** 2 - (np.abs(fft(v_dc - self.Vdc) / N)[0]) ** 2)
        I_dc_eff = np.sqrt(1 / self.Tel / K * np.sum(i_dc ** 2 * dt))
        I_dc_v1_eff = np.abs(fft(i_dc) / N)[0]
        I_dc_thd = np.sqrt(I_dc_eff ** 2 - I_dc_v1_eff ** 2)

        # ==============================================================================
        # Post-Processing
        # ==============================================================================
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
        outDc['V_dc_eff'] = V_dc_eff
        outDc['V_dc_v1_eff'] = V_dc_v1_eff
        outDc['V_dc_thd'] = V_dc_thd
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
    def calcFreq(self, s, xs, i_a, v_a, v_a0, i_dc, v_dc):
        # ==============================================================================
        # Description
        # ==============================================================================
        """
        This function calculates the frequency domain results based on the time domain
        waveforms of the B2 bridge.

        Input:
        1) t:       input time vector (sec)
        2) i_a:     load current (A)
        3) v_a:     load voltage (V)
        4) i_dc:    dc current (A)
        5) v_dc:    dc voltage (V)

        Output:
        1) outAc:   outputs distortion ac side
        2) outDc:   outputs distortion dc side
        """

        # ==============================================================================
        # Initialisation
        # ==============================================================================
        N = int(len(s))
        freqAc = {}
        freqDc = {}
        freqSw = {}

        # ==============================================================================
        # Calculation
        # ==============================================================================
        # ------------------------------------------
        # AC Side
        # ------------------------------------------
        # Sequence
        Y = np.abs(fft(s) / N)[0:int(N / 2)]
        Y[1:-2] = 2 * Y[1:-2]
        freqSw['Sa'] = Y

        # Sampled Reference
        Y = np.abs(fft(xs) / N)[0:int(N / 2)]
        Y[1:-2] = 2 * Y[1:-2]
        freqSw['Xas'] = Y

        # Current
        Y = np.abs(fft(i_a) / N)[0:int(N / 2)]
        Y[1:-2] = 2 * Y[1:-2]
        freqAc['I_a'] = Y

        # Line Voltage
        Y = np.abs(fft(v_a) / N)[0:int(N / 2)]
        Y[1:-2] = 2 * Y[1:-2]
        freqAc['V_a'] = Y

        # Line-Neutral Voltage
        Y = np.abs(fft(v_a0) / N)[0:int(N / 2)]
        Y[1:-2] = 2 * Y[1:-2]
        freqAc['V_a0'] = Y

        # ------------------------------------------
        # DC Side
        # ------------------------------------------
        # Current
        Y = np.abs(fft(i_dc) / N)[0:int(N / 2)]
        Y[1:-2] = 2 * Y[1:-2]
        freqDc['I_dc'] = Y

        # Voltage
        Y = np.abs(fft(v_dc) / N)[0:int(N / 2)]
        Y[1:-2] = 2 * Y[1:-2]
        freqDc['V_dc'] = Y

        # ==============================================================================
        # Return
        # ==============================================================================
        return [freqSw, freqAc, freqDc]

#######################################################################################################################
# References
#######################################################################################################################
