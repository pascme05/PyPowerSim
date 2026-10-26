#######################################################################################################################
#######################################################################################################################
# Title:        Power Electronics Simulation Toolkit (PyPowerSim)
# Topic:        Power Electronics
# File:         tutorials.py
# Date:         26.10.2026
# Author:       Dr. Pascal A. Schirmer
# Version:      V.1.0
# Copyright:    Pascal Schirmer
#######################################################################################################################
#######################################################################################################################

#######################################################################################################################
# Function Description
#######################################################################################################################
"""
Tutorials of the GUI: guided scenarios per operating mode. Every tutorial is a config (config/tutorials/<id>.json)
made of a base config (config/<base>.json: default for the inverter, default_DCDC for the DAB) and the changes listed
below; it uses the device workbooks of the repository. The 'look' texts were checked with runs of the toolkit.
"""

#######################################################################################################################
# Import libs
#######################################################################################################################
from gui import adapter as A
from gui.core.config import CONFIG_DIR, from_json, to_json, default_config, read_config

import copy
import os

TUT_DIR = os.path.join(CONFIG_DIR, 'tutorials')
RUN = {'save': 0, 'plotGen': 0}                                       # run settings of all tutorials

#######################################################################################################################
# Tutorials
#######################################################################################################################
TUTORIALS = [
    # ==================================================================================================================
    # Mode 0: sweep
    # ==================================================================================================================
    {'id': 'tut01_sweep', 'mode': 0, 'level': 'Beginner', 'minutes': 5, 'base': 'default',
     'title': 'Modulation index sweep of a three-phase inverter',
     'story': "A three-phase bridge (B6) at Vdc = 600 V feeds an RL load (R = 5 Ω, L = 5 mH, 50 Hz) with space vector "
              "PWM at fs = 1050 Hz (pulse number 21). How do the output voltage, the current and their distortion "
              "change with the modulation index, and how well does the analytical solution of the toolkit match the "
              "simulation?",
     'goals': ["run a sweep and read RMS, fundamental and distortion over Mi",
               "compare the numerical simulation with the analytical solution",
               "relate the results to the PWM theory (linear range, side bands)"],
     'changes': {'Exp.Name': 'tut01'},
     'steps': ["Open **0 · Sweep › Theory › PWM explorer** and move the modulation index: the fundamental of v_a0 "
               "follows Mi·Vdc/2.",
               "**Parameters**: the sweep covers the modulation index up to Mi = 1 (Dat › stat › Mi) in 20 points.",
               "**Run** the sweep.",
               "**Results**: select V_a_v1_eff and V_a_eff — fundamental and RMS of the phase voltage over Mi, "
               "solid = simulation, dashed = analytical.",
               "Select I_a_thd and V_a_thd: the RMS value of the harmonics (distortion) in A and V."],
     'look': "The run takes about 35 s (20 operating points). The fundamental of the phase voltage grows linearly "
             "with Mi: 111.5 V RMS at Mi = 0.53 and 211 V at Mi = 1 (analytical 111.6 V / 211.9 V = Mi·Vdc/(2√2)). "
             "The RMS value of the PWM voltage is much larger (256 V at Mi = 1) because the harmonics carry "
             "V_a_thd ≈ 146–149 V over the whole range. The load filters them: at Mi = 1 the current has a "
             "fundamental of 40.3 A and only I_a_thd = 2.5 A of harmonics (≈ 6 % of the fundamental). Simulation and "
             "analytical solution agree for fundamentals and distortion; the analytical RMS value of the current "
             "(I_a_eff, dashed) is higher than the simulated one (49.1 A vs 40.3 A at Mi = 1) — it does not equal "
             "√(I_1² + I_thd²) of the analytical columns.",
     'exercise': ["Set **Par › PWM › fs** to 2100 Hz (q = 42) and repeat: the current distortion I_a_thd should "
                  "drop (the inductance filters higher frequencies better), the voltage distortion hardly changes "
                  "(the PWM voltage still switches between the same levels).",
                  "Set **Top › L** to 1 mH: the current distortion grows (less filtering), the voltages stay the same."]},
    # ==================================================================================================================
    # Mode 1: steady state
    # ==================================================================================================================
    {'id': 'tut02_steady', 'mode': 1, 'level': 'Beginner', 'minutes': 10, 'base': 'default',
     'title': 'Losses and temperatures of an IGBT inverter',
     'story': "The same B6 inverter uses the 1200 V / 75 A IGBT IKQ75N120CS6 and an electrolytic DC-link capacitor. "
              "At full modulation (Mi = 1) the inverter delivers about 8 kW into the RL load. How large are the "
              "losses of a switch, how are they split between transistor and diode and conduction and switching, "
              "and how hot do the junctions get?",
     'goals': ["read losses, efficiency and temperatures of the steady state",
               "distinguish conduction and switching losses of transistor and diode",
               "relate the temperature ripple to the thermal network of the device"],
     'changes': {'Exp.Name': 'tut02'},
     'steps': ["**Inputs**: look at the IGBT IKQ75N120CS6 — *Characteristics* (v_CE, E_on, E_off over the current and "
               "Tj) and *Thermal network* (Zth).",
               "**Theory › Losses of a switch**: the analytical loss estimate for the operating point.",
               "**Run** the steady state (it iterates until the temperatures change less than the tolerance).",
               "**Results › Losses**: losses per switch (T: transistor, D: diode) split into conduction and switching.",
               "**Results › Temperatures**: junction temperatures of the transistors and diodes over the last period."],
     'look': "The run takes about 70 s. Output power 8.16 kW, total losses 193 W, efficiency 97.6 %. All six switches "
             "have about 32.2 W: transistor conduction 25.0 W, transistor switching 2.7 W, diode conduction 3.6 W "
             "and diode recovery 0.77 W — at fs = 1050 Hz the conduction losses dominate (the DC-link capacitor "
             "shows ≈ 0 W). The transistor junctions reach 57.6 °C on average and swing between 53.5 and 63.3 °C "
             "within the 50 Hz period (≈ 10 K ripple, case at 25 °C); the diodes stay at 30–34 °C. THD of i_a 6.3 %.",
     'exercise': ["Double the switching frequency (**Par › PWM › fs** = 2100 Hz): the switching losses roughly "
                  "double, the conduction losses stay — compare the shares in Results › Losses.",
                  "Set **Dat › stat › Tc** (case/coolant temperature) to 80 °C: the device characteristics are "
                  "evaluated at the higher junction temperature — compare losses and temperatures."]},
    {'id': 'tut03_dpwm', 'mode': 1, 'level': 'Intermediate', 'minutes': 10, 'base': 'default',
     'title': 'Discontinuous PWM: fewer switching losses',
     'story': "With discontinuous PWM (DPWM1) each leg is clamped to the positive or negative DC rail for 120° of the "
              "fundamental period, where it does not switch. How much switching loss does this save compared with "
              "the space vector PWM of the previous tutorial, and what does it cost in current quality?",
     'goals': ["switch the PWM to a carrier-based method with a zero sequence",
               "compare losses and THD with tutorial 'Losses and temperatures of an IGBT inverter'",
               "see the clamping intervals in the switching functions"],
     'changes': {'Exp.Name': 'tut03', 'Par.PWM.type': 'CB', 'Par.PWM.zero': 'DPWM1'},
     'steps': ["**Theory › PWM explorer** (mode 0): choose PWM CB and zero sequence DPWM1 — the modulation function "
               "x_a is clamped to ±1 around the peaks of the voltage.",
               "**Parameters**: Par › PWM › type = CB, zero = DPWM1, all other settings as in tutorial 2.",
               "**Run** the steady state.",
               "**Results › Losses** and the key figures: compare with tutorial 2.",
               "**Results › Waveforms › Switching**: the modulation function xA sits at ±1 and the switching "
               "function sA stays constant during the clamping intervals."],
     'look': "The run takes about 70 s. xA is clamped to ±1 for one third of the period; leg a switches 30 instead "
             "of 42 times per period. The transistor switching losses drop from 2.7 W to about 2.1 W per switch, the "
             "conduction losses rise slightly (25.4 W instead of 25.0 W), in total 190 W instead of 193 W (efficiency "
             "97.69 % instead of 97.63 %). The current quality gets worse: THD of i_a 8.0 % instead of 6.3 %. The "
             "maximum junction temperature hardly changes (63–64 °C). At 1050 Hz the switching losses are only "
             "≈ 10 % of the losses, so discontinuous PWM pays off mainly where switching losses dominate (higher fs).",
     'exercise': ["Try DPWM0 and DPWM2: the clamping interval shifts relative to the current peak — compare the "
                  "switching losses.",
                  "Run DPWM1 at fs = 1575 Hz: about the same number of switching events as SVPWM at 1050 Hz — "
                  "compare losses and THD with tutorial 2."]},
    {'id': 'tut04_dab', 'mode': 1, 'level': 'Intermediate', 'minutes': 10, 'base': 'default_DCDC',
     'title': 'Dual active bridge: power transfer and losses',
     'story': "A dual active bridge converts 800 V to 400 V (20 kW, 100 kHz) with SiC MOSFET modules and a "
              "transformer with N = 1.78. The output voltage is controlled (Exp › output = V); the toolkit computes "
              "the phase shift from the DAB power equation. Where are the losses, and do the bridges switch at zero "
              "voltage?",
     'goals': ["set up a DC-DC (DAB) simulation and its input files", "relate phase shift and power transfer",
               "read the losses of both bridges and the transformer"],
     'changes': {'Exp.Name': 'tut04'},
     'steps': ["**Setup**: topology DAB — the DC-DC converter input files (primary/secondary switches, capacitors, "
               "transformer) are used, the inverter files not.",
               "**Theory › DAB power transfer** (mode 1): power over the phase shift, current waveform and ZVS.",
               "**Inputs**: transformer trafoDAB — turns ratio and leakage inductance.",
               "**Run** the steady state (about a minute: fsim = 100 MHz).",
               "**Results**: key figures, **Losses** per switch and transformer, **Waveforms** of the transformer "
               "current."],
     'look': "The run takes about 45 s. The toolkit sets the phase shift to 0.878 rad (50.3°) for Vo = 400 V at "
             "R = 8 Ω (Po = 20 kW, as in the theory explorer). Output power 19.9 kW, total losses 532 W, efficiency "
             "97.3 %. Primary switches S1–S4 about 61 W each (switching 37.7 W, diode recovery 9.4 W, conduction "
             "≈ 14 W), secondary switches S5–S8 about 41–42 W, transformer 117 W (48 W and 56 W in the windings, "
             "13 W core), capacitors 1.2 W and 2.5 W. The switching losses dominate although the explorer shows "
             "ZVS for both bridges: the loss model of the toolkit assumes hard switching (Par › PWM › sw = 0). The "
             "primary transistors reach 152 °C (diodes 125 °C), the secondary transistors 110 °C; at 100 kHz the "
             "temperatures show no ripple.",
     'exercise': ["In **Theory › DAB power transfer** set V2 = 300 V and adjust φ: watch the ZVS indicators and the "
                  "RMS current — then run the simulation with **Dat › stat › Vo** = 300 V.",
                  "Set the switching frequency to 50 kHz in the explorer: the same power needs a smaller phase "
                  "shift (P ∝ 1/fs at constant φ)."]},
    # ==================================================================================================================
    # Mode 2: transient
    # ==================================================================================================================
    {'id': 'tut05_transient', 'mode': 2, 'level': 'Intermediate', 'minutes': 10, 'base': 'default',
     'title': 'Heating of the inverter after switch-on',
     'story': "The inverter of tutorial 2 starts at 25 °C. How fast do the junctions heat up, and is the simulated "
              "time of 0.5 s enough to reach the steady state?",
     'goals': ["run a transient simulation", "relate the heating to the time constants of the Foster network",
               "compare with the steady state"],
     'changes': {'Exp.Name': 'tut05'},
     'steps': ["**Theory › Heating of the devices**: Zth step response and the settling of each branch at tmax.",
               "**Parameters**: Dat › trans › tmax = 0.5 s, initial temperatures 25 °C.",
               "**Run** the transient simulation.",
               "**Results › Temperatures** (evaluation window *Full simulation*): heating of transistors and diodes.",
               "Compare the end temperatures with tutorial 2 (steady state)."],
     'look': "The run takes about 2.5 min. The transistor junctions heat from 25 °C to a maximum of 63–65 °C within "
             "the 0.5 s, the diodes to 33–34 °C; in the last period the transistors swing between 53 and 65 °C "
             "(mean 57.4–59.1 °C), practically the steady state of tutorial 2 (57.6 °C mean, 63.3 °C maximum). "
             "0.5 s is enough because the dominant branch case → ambient (1 K/W, τ = 0.1 s) has settled; the slowest "
             "junction → case branch (τ ≈ 3 s) has only 0.0013 K/W. Losses 193 W over the whole run, efficiency "
             "97.6 %, THD of i_a 6.3 %. The switches S2, S4, S6 show about 1 W more switching losses than S1, S3, S5 "
             "in this run.",
     'exercise': ["Shorten tmax to 0.1 s (one τ of the case → ambient branch): the junctions reach only about 70 % "
                  "of the temperature rise (fast junction → case branches settled, case → ambient at 63 %) — check "
                  "with **Theory › Heating of the devices** first.",
                  "Halve the fundamental frequency (**Top › fel** = 25 Hz): the temperature ripple per period grows "
                  "(see **Theory › Temperature ripple** of mode 1)."]},
    # ==================================================================================================================
    # Mode 3: closed loop
    # ==================================================================================================================
    {'id': 'tut06_hysteresis', 'mode': 3, 'level': 'Advanced', 'minutes': 10, 'base': 'default',
     'title': 'Hysteresis current control',
     'story': "Instead of a PWM with fixed switching frequency, a hysteresis controller switches a leg whenever the "
              "current leaves a band around the reference. The reference steps from 50 % to 100 % at half of the "
              "simulated time. How well is the reference tracked and how often do the switches switch?",
     'goals': ["run the closed loop with a short simulated time",
               "see the effect of the band width and the controller frequency fc",
               "read the variable switching frequency in the spectrum"],
     'changes': {'Exp.Name': 'tut06', 'Exp.plot': 1, 'Par.Cont.type': 'HY', 'Par.Cont.hys': 5,
                 'Dat.trans.tmax': 0.04},
     'steps': ["**Theory › Hysteresis current control**: change band and fc and watch the switching frequency.",
               "**Parameters**: Par › Cont › type = HY, band 5 %, fc = 20 kHz; Dat › trans › tmax = 0.04 s (two "
               "periods — the closed loop simulates tmax·fc = 800 control steps one after another).",
               "**Run** the closed loop.",
               "**Results › Waveforms** (window *Full simulation*): i_a and the reference i_ref_a with the step.",
               "**Results › Spectrum**: no distinct switching side bands — the switching frequency varies."],
     'look': "The run takes about 2 min (800 control steps). The reference steps from 28.6 A to 57.2 A peak at "
             "20 ms; the current follows with an RMS error of 1.9 A before and 2.8 A after the step (maximum 9.7 A). "
             "The leg switches on average at 3.2 kHz: 4.6 kHz at half reference, 1.7 kHz at full reference (the band "
             "scales with the reference, the RL load limits di/dt). In the last period the output power is 7.94 kW "
             "(I_rms = 39.9 A), the losses 132 W (efficiency 98.3 %) and the THD of i_a 6.8 % — similar to SVPWM in "
             "tutorial 2 at a lower switching frequency, but the losses differ from switch to switch (11.5 … 29.3 W) "
             "because the switching pattern is not symmetric. The closed loop has no thermal model (no "
             "temperatures).",
     'exercise': ["Set the band to 1 %: the current follows more closely, but fc limits the switching (the current "
                  "overshoots by up to Vdc/(2L·fc)).",
                  "Increase fc to 50 kHz (longer run) and compare the tracking error."]},
]


#######################################################################################################################
# Functions
#######################################################################################################################
def merge_changes(cfg, t):
    """Config of a tutorial: base config + run settings + changes 'Section.cat.var' (or 'Exp.var')."""
    cfg = copy.deepcopy(cfg)
    cfg[A.RUN_SECTION].update(RUN)
    for k, v in t['changes'].items():
        *path, last = k.split('.')
        d = cfg
        for p in path:
            d = d.setdefault(p, {})
        d[last] = copy.deepcopy(v)
    cfg[A.RUN_SECTION][A.MODE_VAR] = t['mode']
    return cfg


def base(t):
    return default_config() if t.get('base', 'default') == 'default' else read_config(t['base'])


def config_file(tid):
    return os.path.join(TUT_DIR, f"{tid}.json")


def load(tid):
    """Config of a tutorial (config/tutorials/<id>.json, generated from base and changes if the file is missing)."""
    f = config_file(tid)
    if os.path.isfile(f):
        with open(f, 'r', encoding='utf-8') as fh:
            return from_json(fh.read())
    return merge_changes(base(get(tid)), get(tid))


def get(tid):
    return next(t for t in TUTORIALS if t['id'] == tid)


def write_all():
    """Writes config/tutorials/<id>.json for all tutorials."""
    os.makedirs(TUT_DIR, exist_ok=True)
    for t in TUTORIALS:
        with open(config_file(t['id']), 'w', encoding='utf-8') as fh:
            fh.write(to_json(merge_changes(base(t), t)))


if __name__ == '__main__':
    write_all()
