#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         classDAB
# Date:         04.02.2026
# Author:       Dr. Pascal A. Schirmer
# Version:      V.1.0
# Copyright:    Pascal Schirmer
#######################################################################################################################
#######################################################################################################################

#######################################################################################################################
# Function Description
#######################################################################################################################
"""
This class initializes an object of a DAB (dual-active-bridge) converter.
Inputs:     1) fel:     electrical frequency at the output (Hz)
            2) fs:      switching frequency of the converter cell (Hz)
            3) fc:      controller frequency (Hz)
            4) fsim:    simulation frequency of the converter cell (Hz)
            5) td:      dead-time of the switching devices (sec)
            6) tmin:    minimum on-time of a pulse (sec)
            7) cyc:     number of cycles till convergence
            8) W:       number of points evaluated for distortion analysis
            9) Mi:      modulation index (p.u.)
            10) phi:    phase shift for DAB (deg)
            11) N:      transformer ratio (Np/Ns)
            12) Lk:     leakage inductance (power transfer inductance) of the DAB (H)
            13) Vdc:    dc link voltage at the input of the converter cell (V)
            14) Tc_st:  case temperature for stationary analysis
            15) Tj_st:  core temperature for stationary analysis
            16) Tc_tr:  case temperature for transient analysis
            17) Tj_tr:  core temperature for transient analysis
"""

#######################################################################################################################
# Import libs
#######################################################################################################################
# ==============================================================================
# Internal
# ==============================================================================
from src.general.helpFnc import deadTime
from src.cont.conHys import conHys

# ==============================================================================
# External
# ==============================================================================
import numpy as np
import pandas as pd
from scipy import signal
import copy


#######################################################################################################################
# Class
#######################################################################################################################
class classDAB:
    ###################################################################################################################
    # Constructor
    ###################################################################################################################
    def __init__(self, fel, fs, fc, fsim, td, tmin, cyc, W, Mi, phi, Ntr, Lk, Vdc, Tc_st, Tj_st, Tc_tr, Tj_tr):
        self.fel = fel
        self.fs = fs
        self.fc = fc
        self.fsim = fsim
        self.td = td
        self.tmin = tmin
        self.K = int(cyc)
        self.W = int(W)
        self.Mi = Mi
        self.phi = phi
        self.Ntr = Ntr
        self.Lk = Lk
        self.Vdc = Vdc
        self.Tc_st = Tc_st
        self.Tj_st = Tj_st
        self.Tc_tr = Tc_tr
        self.Tj_tr = Tj_tr
        self.Nsim = int(np.ceil(fsim / fel))
        self.Npwm = int(np.ceil(fs / fel))
        self.q = int(fs / fel)
        self.N = int(fsim / fel)
        self.Ts = 1 / fs
        self.Tel = 1 / fel
        self.Mi_max = 1
        self.id1 = ['A', 'B', 'C', 'D']
        self.id2 = ['S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8']
        self.id3 = ['A', 'A', 'B', 'B', 'C', 'C', 'D', 'D']
        self.id4 = ['i_ac_pri', 'i_ac_pri', 'i_ac_pri', 'i_ac_pri',
                    'i_ac_sec', 'i_ac_sec', 'i_ac_sec', 'i_ac_sec']
        self.id5 = ['HS', 'LS', 'HS', 'LS', 'HS', 'LS', 'HS', 'LS']
        self.id6 = ['T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8']
        self.id7 = ['D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8']
        self.id8 = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8']
        self.id9 = [+1, +1, -1, -1, +1, +1, -1, -1]
        self.name = 'DAB'

    ###################################################################################################################
    # Init Data
    ###################################################################################################################
    def initData(self):
        # ==============================================================================
        # Description
        # ==============================================================================
        """
        This function initializes the data structure of the DAB bridge.

        Input:

        Output:
        1) data:        initialised data structure
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
        for sid in self.id2:
            data['elec']['sw'][sid] = pd.DataFrame(columns=['i_T', 'v_T', 'i_D', 'v_D'])

        # ------------------------------------------
        # Losses
        # ------------------------------------------
        data['loss']['sw'] = {}
        for sid in self.id2:
            data['loss']['sw'][sid] = pd.DataFrame(columns=['p_T_c', 'p_T_s', 'p_T', 'p_D_c', 'p_D_s', 'p_D', 'p_L'])

        # ------------------------------------------
        # Thermal
        # ------------------------------------------
        data['ther']['sw'] = pd.DataFrame(columns=self.id6 + self.id7 + self.id8)

        # ==============================================================================
        # Capacitor
        # ==============================================================================
        # ------------------------------------------
        # Electric
        # ------------------------------------------
        data['elec']['cap'] = {}
        data['elec']['cap']['C1'] = pd.DataFrame(columns=['i_c', 'v_c'])
        data['elec']['cap']['C2'] = pd.DataFrame(columns=['i_c', 'v_c'])

        # ------------------------------------------
        # Losses
        # ------------------------------------------
        data['loss']['cap'] = {}
        data['loss']['cap']['C1'] = pd.DataFrame(columns=['p_L'])
        data['loss']['cap']['C2'] = pd.DataFrame(columns=['p_L'])

        # ------------------------------------------
        # Thermal
        # ------------------------------------------
        data['ther']['cap'] = pd.DataFrame(columns=['C1', 'C2'])

        # ==============================================================================
        # Transformer
        # ==============================================================================
        # ------------------------------------------
        # Electric
        # ------------------------------------------
        data['elec']['tra'] = {}
        data['elec']['tra']['T1'] = pd.DataFrame(columns=['i_p', 'i_s', 'i_m', 'v_p', 'v_s', 'Phi', 'B'])

        # ------------------------------------------
        # Losses
        # ------------------------------------------
        data['loss']['tra'] = {}
        data['loss']['tra']['T1'] = pd.DataFrame(columns=['p_L_pri', 'p_L_sec', 'p_L_core', 'p_L'])

        # ------------------------------------------
        # Thermal
        # ------------------------------------------
        data['ther']['tra'] = {}
        data['ther']['tra']['T1'] = pd.DataFrame(columns=['Pri', 'Sec', 'Core'])

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
        This function initializes the output files of the B6 bridge

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
        distAc = {'Pri': {}, 'Sec': {}}
        distDc = {'Pri': {}, 'Sec': {}}
        timeElec = {'sw': {}, 'cap': {}, 'tra': {}}
        timeLoss = {'sw': {}, 'cap': {}, 'tra': {}}
        timeTher = {'sw': {}, 'cap': {}, 'tra': {}}

        # ==============================================================================
        # Calculation
        # ==============================================================================
        # ------------------------------------------
        # Time
        # ------------------------------------------
        timeElec['cap']['C1'] = pd.DataFrame(columns=['i_c', 'v_c'])
        timeElec['cap']['C2'] = pd.DataFrame(columns=['i_c', 'v_c'])
        timeSw = pd.DataFrame(columns=['t', 'v_a_ref', 'v_b_ref', 'e', 'xAs', 'xBs', 'xAsh', 'xBsh',
                                       'sA', 'sB', 'sC', 'sD', 'cA', 'cB', 'c'])
        
        # ------------------------------------------
        # Frequency
        # ------------------------------------------
        freqSw = pd.DataFrame(columns=['Sa', 'Sb', 'Xas', 'Xbs'])
        freqAc = pd.DataFrame(columns=['I_ab', 'V_ab'])
        freqDc = pd.DataFrame(columns=['I_dc', 'I_d_p', 'I_d_m', 'V_dc'])
        
        # ------------------------------------------
        # Distortion
        # ------------------------------------------
        # AC
        distAc['Pri']['num'] = pd.DataFrame(data=np.ones((self.W, 6)),
                                            columns=['V_a_eff', 'V_a_v1_eff', 'V_a_thd', 'I_a_eff', 'I_a_v1_eff', 'I_a_thd'])
        distAc['Sec']['num'] = pd.DataFrame(data=np.ones((self.W, 6)),
                                            columns=['V_a_eff', 'V_a_v1_eff', 'V_a_thd', 'I_a_eff', 'I_a_v1_eff', 'I_a_thd'])
        distAc['Pri']['ana'] = pd.DataFrame(data=np.ones((self.W, 6)),
                                            columns=['V_a_eff', 'V_a_v1_eff', 'V_a_thd', 'I_a_eff', 'I_a_v1_eff', 'I_a_thd'])
        distAc['Sec']['ana'] = pd.DataFrame(data=np.ones((self.W, 6)),
                                            columns=['V_a_eff', 'V_a_v1_eff', 'V_a_thd', 'I_a_eff', 'I_a_v1_eff', 'I_a_thd'])
        
        # DC
        distDc['Pri']['num'] = pd.DataFrame(data=np.ones((self.W, 6)),
                                            columns=['V_dc_eff', 'V_dc_v1_eff', 'V_dc_thd', 'I_dc_eff', 'I_dc_v1_eff', 'I_dc_thd'])
        distDc['Sec']['num'] = pd.DataFrame(data=np.ones((self.W, 6)),
                                            columns=['V_dc_eff', 'V_dc_v1_eff', 'V_dc_thd', 'I_dc_eff', 'I_dc_v1_eff', 'I_dc_thd'])
        distDc['Pri']['ana'] = pd.DataFrame(data=np.ones((self.W, 6)),
                                            columns=['V_dc_eff', 'V_dc_v1_eff', 'V_dc_thd', 'I_dc_eff', 'I_dc_v1_eff', 'I_dc_thd'])
        distDc['Sec']['ana'] = pd.DataFrame(data=np.ones((self.W, 6)),
                                            columns=['V_dc_eff', 'V_dc_v1_eff', 'V_dc_thd', 'I_dc_eff', 'I_dc_v1_eff', 'I_dc_thd'])

        # ==============================================================================
        # Return
        # ==============================================================================
        return [timeSw, timeElec, timeLoss, timeTher, freqSw, freqDc, freqAc, distAc, distDc]

    ###################################################################################################################
    # Init Output
    ###################################################################################################################
    def out(self, timeElec, timeLoss, timeTher, timeAc, timeDc, freqSw, freqAc, freqDc, distAc, distDc, t_ref, v_ref,
            e_ref, s, c, xs, xsh, x, xN0, M_i, t0, t1, Nel):
        """
        This function summarizes the output structure of the DAB bridge.
        """

        time = {}
        freq = {}
        sweep = {}
        [timeSw, _, _, _, _, _, _, _, _] = self.initOut()

        # Align lengths across all vectors to avoid index mismatch
        def _len(arr):
            return len(arr[t0:t1])

        n = min(
            _len(t_ref),
            _len(v_ref['A']), _len(v_ref['B']),
            _len(e_ref['A']),
            _len(s['A']), _len(s['B']), _len(s['C']), _len(s['D']),
            _len(c['A']), _len(c['B']),
            _len(xs['A']), _len(xs['B']),
            _len(xsh['A']), _len(xsh['B']),
            _len(x['A']), _len(x['B']),
            _len(xN0['A']), _len(xN0['B'])
        )

        t0n = t0
        t1n = t0 + n

        timeSw['t'] = t_ref[t0n:t1n]
        timeSw['v_a_ref'] = v_ref['A'][t0n:t1n]
        timeSw['v_b_ref'] = v_ref['B'][t0n:t1n]
        timeSw['e_a'] = e_ref['A'][t0n:t1n]
        timeSw['sA'] = s['A'][t0n:t1n]
        timeSw['sB'] = s['B'][t0n:t1n]
        timeSw['sC'] = s['C'][t0n:t1n]
        timeSw['sD'] = s['D'][t0n:t1n]
        timeSw['c'] = c['A'][t0n:t1n]
        timeSw['cA'] = c['A'][t0n:t1n]
        timeSw['cB'] = c['B'][t0n:t1n]
        timeSw['xAs'] = xs['A'][t0n:t1n]
        timeSw['xBs'] = xs['B'][t0n:t1n]
        timeSw['xAsh'] = xsh['A'][t0n:t1n]
        timeSw['xBsh'] = xsh['B'][t0n:t1n]
        timeSw['xA'] = x['A'][t0n:t1n]
        timeSw['xB'] = x['B'][t0n:t1n]
        timeSw['n0'] = xN0['A'][t0n:t1n]
        timeSw['n0A'] = xN0['A'][t0n:t1n]
        timeSw['n0B'] = xN0['B'][t0n:t1n]

        if not timeLoss:
            time['t'] = timeSw['t']
        else:
            time['t'] = np.linspace(0, self.Tel * Nel, int(len(timeLoss['sw'][self.id2[0]]['p_T'])))

        time['Sw'] = timeSw
        time['Ac'] = timeAc
        time['Dc'] = timeDc
        time['Elec'] = timeElec
        time['Loss'] = timeLoss
        time['Ther'] = timeTher

        freq['Sw'] = freqSw
        freq['Ac'] = freqAc
        freq['Dc'] = freqDc

        sweep['Ac'] = distAc
        sweep['Dc'] = distDc
        sweep['Mi'] = M_i

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
    # PWM / Phase-Shift Control
    ###################################################################################################################
    def calcPWM(self, _1, t_ref, Mi, _2):
        # ==============================================================================
        # Description
        # ==============================================================================
        """
        This function calculates the switching times assuming any provided PWM method.
        Currently, only SPS is implemented. The modulation index is always 0.5.

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
        if len(t_ref) == 0:
            t_ref = np.linspace(0, self.K / self.fel, self.K * self.N + 1)
        phase = 2 * np.pi * self.fs * t_ref

        # ==============================================================================
        # Calculation
        # ==============================================================================
        s_p = np.where(np.sin(phase) >= 0, 1, -1)
        s_s = np.where(np.sin(phase - self.phi) >= 0, 1, -1)
        sA = self.calcDead(s_p.copy(), t_ref, Mi)
        sB = -sA
        sC = self.calcDead(s_s.copy(), t_ref, Mi)
        sD = -sC

        # ==============================================================================
        # Output
        # ==============================================================================
        s = {'A': sA, 'B': sB, 'C': sC, 'D': sD}
        c = {'A': sA, 'B': sB}
        xs = {'A': sA, 'B': sB}
        xsh = {'A': sA, 'B': sB}
        x = {'A': sA, 'B': sB}
        xN0 = {'A': np.zeros(np.size(sA)), 'B': np.zeros(np.size(sB))}

        # ==============================================================================
        # Return
        # ==============================================================================
        return [xs, xsh, s, c, x, xN0]

    ###################################################################################################################
    # Closed Loop Control
    ###################################################################################################################
    def calcCON(self, i_ref, i_act, s_act, t_con, scale, setup):
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
        4) t_con:   time instance of the control action (sample)
        5) scale:   scaling value to create step response
        5) setup:   variable including all parameters

        Output:
        1) s:       switching instances (sec)
        2) Mi:      updated modulation index
        3) err:     error between reference and actual signal
        """

        # ==============================================================================
        # Init
        # ==============================================================================
        s_out = {}

        # Helper: safely pick a signal value from a dict at index t_con
        def _pick_signal(sig_dict, keys, idx):
            for key in keys:
                if key in sig_dict:
                    sig = sig_dict[key]
                    break
            else:
                return 0.0
            try:
                sig = np.asarray(sig)
                if sig.size == 0:
                    return 0.0
                return float(sig[idx])
            except Exception:
                try:
                    return float(sig)
                except Exception:
                    return 0.0

        # Helper: get full signal for tolerance calculation
        def _pick_signal_arr(sig_dict, keys):
            for key in keys:
                if key in sig_dict:
                    sig = sig_dict[key]
                    break
            else:
                return np.array([0.0])
            try:
                return np.asarray(sig)
            except Exception:
                try:
                    return np.array([float(sig)])
                except Exception:
                    return np.array([0.0])

        # Reference amplitude for hysteresis band
        ref_for_tol = _pick_signal_arr(i_ref, ['i_dc_out', 'A', 'i_ac_pri', 'i_ac_sec'])
        tol = np.max(np.abs(ref_for_tol)) * setup['Par']['Cont']['hys'] / 100 * scale

        # ==============================================================================
        # Calculation
        # ==============================================================================
        i_act_val = _pick_signal(i_act, ['i_dc_out', 'i_ac_pri', 'i_ac_sec', 'i_a'], t_con)
        i_ref_val = _pick_signal(i_ref, ['i_dc_out', 'A', 'i_ac_pri', 'i_ac_sec'], t_con) * scale

        # ==============================================================================
        # Calculation
        # ==============================================================================
        if setup['Par']['Cont']['type'] == "HY":
            [s, err] = conHys(i_act_val, i_ref_val, s_act['A'], tol)
        elif setup['Par']['Cont']['type'] == "PI":
            # PI controller to set phase shift (phi) for DAB
            Kp = setup['Par']['Cont'].get('Kp', 0.0)
            Ki = setup['Par']['Cont'].get('Ki', 0.0)
            dt = 1 / self.fc if self.fc != 0 else 0.0

            # Initialize PI integrator if needed
            if not hasattr(self, '_pi_int'):
                self._pi_int = 0.0
            if not hasattr(self, '_pwm_phase'):
                self._pwm_phase = 0.0

            err = i_ref_val - i_act_val
            # Integrate
            self._pi_int = self._pi_int + err * dt
            u_unsat = Kp * err + Ki * self._pi_int

            # Saturation (default: +/- 90 deg)
            phi_lim = np.pi / 2
            sat_min = setup['Par']['Cont'].get('satMin', -phi_lim)
            sat_max = setup['Par']['Cont'].get('satMax', phi_lim)
            try:
                phi_min = max(-phi_lim, float(sat_min))
                phi_max = min(phi_lim, float(sat_max))
            except Exception:
                phi_min, phi_max = -phi_lim, phi_lim

            u_sat = np.clip(u_unsat, phi_min, phi_max)

            # Anti-windup: freeze integrator if saturated and error drives further into saturation
            if (u_unsat != u_sat) and Ki != 0:
                if (u_sat >= phi_max and err > 0) or (u_sat <= phi_min and err < 0):
                    self._pi_int = self._pi_int - err * dt
                    u_unsat = Kp * err + Ki * self._pi_int
                    u_sat = np.clip(u_unsat, phi_min, phi_max)

            phi_cmd = u_sat

            # Build switching sequences for the current control period
            Ncon = len(s_act['A'])
            t_ref = np.arange(Ncon) / self.fsim
            phase = self._pwm_phase + 2 * np.pi * self.fs * t_ref
            s_p = np.where(np.sin(phase) >= 0, 1, -1)
            s_s = np.where(np.sin(phase - phi_cmd) >= 0, 1, -1)

            # Apply dead-time/min pulse
            Mi_loc = setup['Dat']['stat'].get('Mi', 0.5)
            sA = self.calcDead(s_p.copy(), t_ref, Mi_loc)
            sC = self.calcDead(s_s.copy(), t_ref, Mi_loc)

            s_out['A'] = sA
            s_out['B'] = -s_out['A']
            s_out['C'] = sC
            s_out['D'] = -s_out['C']

            # Update phase accumulator for continuity
            self._pwm_phase = (self._pwm_phase + 2 * np.pi * self.fs * (Ncon / self.fsim)) % (2 * np.pi)

            Mi = abs(phi_cmd) / (np.pi / 2) if (np.pi / 2) != 0 else 0.0
            return [s_out, Mi, err]
        else:
            [s, err] = conHys(i_act_val, i_ref_val, s_act['A'], tol)

        # ==============================================================================
        # Post-Processing
        # ==============================================================================
        s_out['A'] = s
        s_out['B'] = -s_out['A']
        s_out['C'] = s_out['A']
        s_out['D'] = -s_out['A']
        Mi = abs(np.mean(s))

        # ==============================================================================
        # Return
        # ==============================================================================
        return [s_out, Mi, err]

    ###################################################################################################################
    # Init controller variables
    ###################################################################################################################
    def initCON(self):
        # ==============================================================================
        # Description
        # ==============================================================================
        """
        This function initializes the relevant variables for calculating closed loop
        control outputs.

        Input:

        Output:
        1) s:       switching instances (sec)
        2) i_act:   actual current vector (A)
        3) swOut:   switching function output
        """

        # ==============================================================================
        # Init
        # ==============================================================================
        Ncon = int(self.fsim / self.fc)

        # ==============================================================================
        # Calculation
        # ==============================================================================
        s = {'A': np.ones(Ncon), 'B': np.ones(Ncon), 'C': np.ones(Ncon), 'D': np.ones(Ncon)}
        i_act = {'i_a': np.zeros(Ncon)}
        outSw = {'A': [], 'B': [], 'C': [], 'D': []}

        # Reset PI controller states for closed-loop runs
        self._pi_int = 0.0
        self._pwm_phase = 0.0

        # ==============================================================================
        # Return
        # ==============================================================================
        return [s, i_act, outSw]

    ###################################################################################################################
    # Init controller variables
    ###################################################################################################################
    def appCON(self, s_i, outSw, iterC):
        # ==============================================================================
        # Description
        # ==============================================================================
        """
        This function appends the calculated controller outputs.

        Input:
        1) s_i:     switching instances (sec)
        2) outSw:   switching sequence output (sec)
        3) iterC:   controller iteration

        Output:
        1) t_con:   time vector for one control iteration (sec)
        2) e_con:   back emf of the controller (V)
        3) outSw:   switching function output
        """

        # ==============================================================================
        # Calculation
        # ==============================================================================
        outSw['A'] = np.append(outSw['A'], s_i['A'])
        outSw['B'] = np.append(outSw['B'], s_i['B'])
        outSw['C'] = np.append(outSw['C'], s_i['C'])
        outSw['D'] = np.append(outSw['D'], s_i['D'])
        e_con = {'A': np.zeros(len(outSw['A']))}
        t_con = np.linspace(0, (iterC + 1) / self.fc, len(outSw['A']))

        # ==============================================================================
        # Return
        # ==============================================================================
        return [outSw, e_con, t_con]

    ###################################################################################################################
    # Temporal Output
    ###################################################################################################################
    def calcTime(self, s, vp, vs, e, t, Mi, mdl, t0, t1, init, para, setup):
        # ==============================================================================
        # Description
        # ==============================================================================
        """
        This function calculates the temporal outputs of the B4 bridge (AC and DC side).

        Input:
        1) s:       switching function
        2) vp:      voltage at the primary bridge (V)
        3) vs:      voltage at the secondary bridge (V)
        2) e:       induced voltage
        3) t:       reference time (sec)
        4) Mi:      modulation index (0 ... 1)
        5) mdl:     transfer functions
        6) t0:      start time (sample)
        7) t1:      end time (sample)
        9) init:    initial conditions for the lsim solver
        10) para:   all parameters used in the simulation
        11) setup:  all setup variables

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
        Lk = para['Tra']['Elec']['con']['Lk']

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
            init = {'inp': [0, 0], 'out': 0, 'dc': 0, 'load': 0}

        # ==============================================================================
        # Calculation
        # ==============================================================================
        # ------------------------------------------
        # AC Voltages
        # ------------------------------------------
        if vp is not None and vs is not None:
            v_ac_pri = vp
            v_ac_sec = vs
            v_ac_sec_ref = self.Ntr * v_ac_sec
            v_ab = v_ac_pri - v_ac_sec_ref
        else:
            v_ac_pri = 0.5 * self.Vdc * (s['A'] - s['B'])
            v_ac_sec = 0.5 * (self.Vdc / self.Ntr) * (s['C'] - s['D']) if self.Ntr != 0 else 0.5 * self.Vdc * (s['C'] - s['D'])
            v_ac_sec_ref = self.Ntr * v_ac_sec
            v_ab = v_ac_pri - v_ac_sec_ref

        # ------------------------------------------
        # AC Currents
        # ------------------------------------------
        _, i_l, _, = signal.lsim(mdl['SS']['AC'], v_ab, t, X0=init['load']*Lk)
        _, i_m, _, = signal.lsim(mdl['SS']['Tra'], v_ac_pri, t, X0=0)
        i_ac_pri = i_l[t0:t1] - np.mean(i_l[t0:t1])
        i_ac_sec = self.Ntr * (i_ac_pri - i_m[t0:t1]) if self.Ntr != 0 else i_ac_pri

        # ------------------------------------------
        # DC-side currents from modulation
        # ------------------------------------------
        i_dc_pri = 0.5 * i_ac_pri * (s['A'][t0:t1] - s['B'][t0:t1])
        i_dc_sec = 0.5 * i_ac_sec * (s['C'][t0:t1] - s['D'][t0:t1])

        # ------------------------------------------
        # Capacitor currents
        # ------------------------------------------
        i_c_pri = i_dc_pri - np.mean(i_dc_pri)
        i_c_sec = i_dc_sec - np.mean(i_dc_sec)

        # ------------------------------------------
        # DC Voltages
        # ------------------------------------------
        # Caps
        _, v_dc_pri_cap, _, = signal.lsim(mdl['SS']['Inp'], i_c_pri, t[t0:t1], X0=init['dc'])
        _, v_dc_sec_cap, _, = signal.lsim(mdl['SS']['Out'], i_c_sec, t[t0:t1], X0=init['dc'])

        # Input Side
        try:
            v_dc_pri = self.Vdc + v_dc_pri_cap - para['CapPri']['Elec']['con']['ESR'] * i_c_pri
        except:
            v_dc_pri = self.Vdc + v_dc_pri_cap

        # Output Side
        try:
            v_dc_sec = np.mean(i_dc_sec) * setup['Top']['R'] + v_dc_sec_cap + para['CapSec']['Elec']['con']['ESR'] * i_c_sec
        except:
            v_dc_sec = np.mean(i_dc_sec) * setup['Top']['R'] + v_dc_sec_cap

        # ------------------------------------------
        # DC Currents
        # ------------------------------------------
        i_dc_in = i_dc_pri - i_c_pri
        i_dc_out = i_dc_sec - i_c_sec

        # ==============================================================================
        # Post-Processing
        # ==============================================================================
        # ------------------------------------------
        # AC Side
        # ------------------------------------------
        # Old for consistency
        outAc['v_ac_pri'] = v_ac_pri[t0:t1]
        outAc['v_ac_sec'] = v_ac_sec[t0:t1]
        outAc['v_ac_sec_ref'] = v_ac_sec_ref[t0:t1]
        outAc['v_L'] = v_ab[t0:t1]
        outAc['i_ac_pri'] = i_ac_pri
        outAc['i_ac_sec'] = i_ac_sec

        # ------------------------------------------
        # DC Side
        # ------------------------------------------
        outDc['v_dc_pri'] = v_dc_pri
        outDc['v_dc_sec'] = v_dc_sec
        outDc['v_dc_pri_cap'] = v_dc_pri_cap
        outDc['v_dc_sec_cap'] = v_dc_sec_cap
        outDc['i_dc_in'] = i_dc_in
        outDc['i_dc_pri'] = i_dc_pri
        outDc['i_dc_sec'] = i_dc_sec
        outDc['i_dc_out'] = i_dc_out
        outDc['i_c_pri'] = i_c_pri
        outDc['i_c_sec'] = i_c_sec

        # ------------------------------------------
        # Init Conditions
        # ------------------------------------------
        outIn['inp'] = 0
        outIn['out'] = 0
        outIn['dc'] = 0
        outIn['load'] = i_ac_pri[-1]

        # ==============================================================================
        # Return
        # ==============================================================================
        return [outAc, outDc, outIn]

    ###################################################################################################################
    # Analytical distortion
    ###################################################################################################################
    def calcDist(self, _1, _2, _3, _4, _5, _6):
        # ==============================================================================
        # Description
        # ==============================================================================
        """
        Placeholder for analytical distortion of the DAB bridge.
        """

        # ==============================================================================
        # Initialisation
        # ==============================================================================
        outAc = {
            'V_a_eff': np.nan, 'V_a_v1_eff': np.nan, 'V_a_thd': np.nan,
            'I_a_eff': np.nan, 'I_a_v1_eff': np.nan, 'I_a_thd': np.nan
        }
        outDc = {
            'V_dc_eff': np.nan, 'V_dc_v1_eff': np.nan, 'V_dc_thd': np.nan,
            'I_dc_eff': np.nan, 'I_dc_v1_eff': np.nan, 'I_dc_thd': np.nan
        }

        # ==============================================================================
        # Return
        # ==============================================================================
        return [outAc, outDc]

    ###################################################################################################################
    # Reference generation
    ###################################################################################################################
    def calcRef(self, _1, _2, _3, t, setup):
        # ==============================================================================
        # Description
        # ==============================================================================
        """
        This function calculates the reference voltage, the reference current functions 
        based on the DAB topology and the given parameters.

        Input:
        1) E:       amplitude of the back emf (V)
        2) phiE:    angle of the back emf (rad)
        3) v_a:     load angle of the output (rad)
        4) t:       given time vector (sec)
        5) setup:   file including all setup parameters

        Output:
        1) v_ref:   reference voltage for a given load scenario (V)
        2) e_ref:   reference back emf for a given load scenario (V)
        3) i_ref:   reference current for a given load scenario (A)
        """

        # ==============================================================================
        # Initialisation
        # ==============================================================================
        v_ref = {}
        e_ref = {}
        i_ref = {}
        if t == [] or len(t) == 0:
            t = np.linspace(0, self.K / self.fel, self.K * self.N + 1)

        # ==============================================================================
        # Calculation
        # ==============================================================================
        # ------------------------------------------
        # Normalized
        # ------------------------------------------
        phase = 2 * np.pi * self.fs * t
        s_p = np.where(np.sin(phase) >= 0, 1, -1)
        s_s = np.where(np.sin(phase - self.phi) >= 0, 1, -1)

        # ------------------------------------------
        # Absolute
        # ------------------------------------------
        v_ref['A'] = 0.5 * self.Vdc * s_p
        v_ref['B'] = 0.5 * (self.Vdc / self.Ntr) * s_s if self.Ntr != 0 else 0.5 * self.Vdc * s_s
        e_ref['A'] = np.zeros(np.size(t))
        e_ref['B'] = np.zeros(np.size(t))
        try:
            Io = setup['Dat']['stat'].get('Io', 0)
            Io = float(np.real(Io))
        except Exception:
            Io = 0.0
        i_ref['A'] = Io * np.ones(np.size(t))
        i_ref['B'] = Io * np.ones(np.size(t))
        i_ref['i_dc_out'] = i_ref['A']

        # ==============================================================================
        # Return
        # ==============================================================================
        return [v_ref, e_ref, i_ref]

#######################################################################################################################
# References
#######################################################################################################################
