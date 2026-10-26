#######################################################################################################################
#######################################################################################################################
# Title:        Power Electronics Simulation Toolkit (PyPowerSim)
# Topic:        Power Electronics
# File:         chart_help.py
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
Explanations of the interactive charts (shown when hovering over the info line above a chart): what the chart shows,
how to read it and what to look for. Keys: name of the view function that calls style.show -> chart key (last part of
the Streamlit key; 'sw' also matches keys like 'res|0|sw|age').
"""

CHART_HELP = {
    # ==================================================================================================================
    # Steady state, transient, closed loop
    # ==================================================================================================================
    '_waveforms': {
        'wave': "**Time signals of the simulation.** Switching: references, carriers and switching functions; AC side: "
                "phase voltages and currents; DC side: DC-link voltage and currents. The evaluation window above "
                "selects the last fundamental period (stationary result) or the whole simulation (start-up, "
                "transient). Zoom with the mouse, click a legend entry to hide a signal.",
    },
    '_loss_view': {
        'bar': "**Mean losses per component over the evaluation window**, split into conduction and switching losses "
               "of transistor and diode (transformer: core and copper losses). Switching losses scale with the "
               "switching frequency, conduction losses with the current squared (MOSFET) or linearly (IGBT, diode).",
    },
    '_temp_view': {
        'temp': "**Temperatures of the components** (T: transistor, D: diode, C: case of a switch). In the stationary "
                "analysis the ripple shows the temperature swing per fundamental period, which drives the thermal "
                "cycling of the devices; in the transient analysis the curves show the heating with the thermal time "
                "constants of the Foster networks.",
    },
    '_spectrum': {
        'spec': "**Amplitude spectrum over the harmonic order** (frequency / fundamental frequency), as in the plots of "
                "the toolkit. Order 1 is the fundamental; the side bands around multiples of fs/fel are the switching "
                "harmonics. The logarithmic scale shows the small low-order harmonics caused by dead time or "
                "saturation.",
    },
    # ==================================================================================================================
    # Sweep
    # ==================================================================================================================
    'pps_sweep': {
        'sweep': "**Sweep of the operating point**: RMS value, fundamental and THD of the voltages and currents over "
                 "the modulation index (inverter) or the phase shift (DAB). Solid lines are the numerical simulation, "
                 "dashed lines the analytical solution of the toolkit — deviations show the effect of dead time, "
                 "sampling and non-ideal devices.",
    },
    # ==================================================================================================================
    # Input data: device workbooks
    # ==================================================================================================================
    '_lut_chart': {
        'Vf': "**On-state voltage of the transistor** v_CE (IGBT) or v_DS (MOSFET) over the current, one curve per "
              "junction temperature. The conduction loss is v·i per sample (src/elec/calcLossSwi.py); a steeper curve "
              "means a higher on-state resistance.",
        'Vfd': "**Forward voltage of the diode** over the diode current, one curve per junction temperature. Used for "
               "the diode conduction losses v_F·i.",
        'Eon': "**Turn-on energy** per switching event over the current, measured at the voltage Vnom of the "
               "workbook; the toolkit scales it linearly with the blocking voltage (Vdc / Vnom) and multiplies with "
               "the number of switching events.",
        'Eoff': "**Turn-off energy** per switching event over the current, measured at Vnom; scaled with Vdc / Vnom.",
        'Erec': "**Reverse-recovery energy of the diode** per switching event over the diode current, measured at "
                "Vnom; scaled with Vdc / Vnom.",
        'Ciss': "**Input capacitance** over the drain-source (collector-emitter) voltage; determines the gate charge "
                "and the switching speed.",
        'Coss': "**Output capacitance** over the blocking voltage; its stored energy is lost at hard turn-on and "
                "enables zero-voltage switching (DAB).",
        'Crss': "**Reverse transfer (Miller) capacitance** over the blocking voltage; sets the voltage slope during "
                "switching (model 'tab' with MOSFETs).",
        'C': "**Capacitance** over the frequency, one curve per temperature.",
        'tan': "**Loss factor tan δ** over the frequency; the ESR follows as tan δ / (2πf·C).",
        'Kr': "**Ripple-current factor** over the frequency: permitted ripple current relative to the rated value.",
        'ESR': "**Equivalent series resistance** over the frequency; the capacitor loss is ESR·I_rms² per harmonic.",
        'Rp': "**Resistance of the primary winding** over the frequency (skin and proximity effect).",
        'Rs': "**Resistance of the secondary winding** over the frequency.",
        'Lm': "**Magnetising inductance** over the flux density (saturation of the core).",
        'Rc1': "**Core-loss resistance** over the flux density.",
        'Rc2': "**Core-loss resistance** over the primary voltage.",
        'Cp': "**Stray capacitance of the primary winding** over the frequency.",
        'Cs': "**Stray capacitance of the secondary winding** over the frequency.",
        'Cps': "**Inter-winding capacitance** primary–secondary over the frequency.",
    },
    'device_thermal': {
        'zth': "**Thermal impedance Zth(t) of the Foster networks** of the thermal sheet: Zth(t) = Σ Rth·(1 − "
               "e^(−t/τ)), τ = Rth·Cth. The temperature rise after a loss step P is P·Zth(t); the end value is the "
               "steady-state thermal resistance ΣRth. Logarithmic time axis to show the fast (chip) and slow (case, "
               "heatsink) branches.",
    },
    # ==================================================================================================================
    # Theory explorers (gui/core/theory.py)
    # ==================================================================================================================
    '_pwm': {
        'pwm': "**One fundamental period computed with the toolkit** (calcRef, calcPWM, calcTime of the topology "
               "class, ideal switches). Top: modulation function, zero sequence and carrier in p.u. of Vdc/2; middle: "
               "leg voltage switching between ±Vdc/2; bottom: phase current filtered by the RL load.",
        'spec': "**Amplitude spectrum of the leg voltage** over the harmonic order (log scale). The fundamental is "
                "order 1; the switching harmonics form side bands around multiples of the pulse number q = fs/fel.",
    },
    '_mi_curve': {
        'mi': "**Fundamental of the leg voltage over the modulation index** (toolkit, one period per point) compared "
              "with linear theory and the six-step limit 4/π. Below the limit of the PWM method the points lie on the "
              "dashed line; above it the gain drops (overmodulation).",
    },
    '_dab_power': {
        'pphi': "**Transferred power over the phase shift** (ideal DAB, single phase shift): parabolic, maximum at "
                "90°. The marker is the selected operating point.",
        'wave': "**Inductor current and bridge voltages over one switching period.** The current changes with slope "
                "(v_1 ± N·v_2)/L_k; its sign at the switching instants decides on zero-voltage switching.",
    },
    '_losses': {
        'bar': "**Mean losses of one switch position** (transistor and anti-parallel diode) over a fundamental period, "
               "from the device characteristics of the workbook and a sinusoidal current.",
        'curves': "**Device characteristics at the selected junction temperature**: on-state voltages (left axis) "
                  "and switching energies scaled to Vdc (right axis) over the current. Shaded: current range of the "
                  "operating point.",
    },
    '_ripple': {
        'tj': "**Junction temperature over three fundamental periods** for a half-wave loss pulse with mean P "
              "(calcTherRC with the Foster network junction → ambient), starting at the mean temperature.",
        'dtf': "**Temperature ripple over the fundamental frequency**: at low frequencies the chip and solder "
               "layers follow the loss pulses (large ΔTj, thermal cycling); at high frequencies the thermal "
               "capacitances average them.",
    },
    '_heating': {
        'heat': "**Heating after a loss step** computed with calcTherRC (solid) and the analytical P·Zth(t) "
                "(dotted, identical). The vertical line marks tmax of the transient mode; the table lists how far "
                "each branch has settled at tmax.",
    },
    '_startup': {
        'i': "**Start-up current of the RL load** for a sinusoidal voltage switched on at angle α: steady sinusoid "
             "plus a DC component decaying with τ = L/R. The DC component vanishes for α = φ.",
    },
    '_hysteresis': {
        'hys': "**Hysteresis current control of one bridge leg** with the controller of the toolkit (conHys), "
               "evaluated once per control period 1/fc. Top: current, reference and band ±δ (reference step from 50 % "
               "to 100 % at half time); bottom: switching state. Sampling lets the current leave the band by up to "
               "Vdc/(2L·fc).",
    },
}


def lookup(scope, key):
    """Explanation of the chart with the Streamlit key `key` drawn by the view function `scope` (or None)."""
    texts = CHART_HELP.get(scope)
    if not texts or not key:
        return None
    parts = str(key).split('|')
    for n in range(len(parts), 0, -1):                     # most specific suffix first: 'exc|curve' -> 'exc'
        for i in range(len(parts) - n, -1, -1):          # from the end: the chart key is the last part
            k = '|'.join(parts[i:i + n])
            if k in texts and i > 0:
                return texts[k]
    return None
