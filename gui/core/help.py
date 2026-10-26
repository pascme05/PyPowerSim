#######################################################################################################################
#######################################################################################################################
# Title:        Power Electronics Simulation Toolkit (PyPowerSim)
# Topic:        Power Electronics
# File:         help.py
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
This function holds the extended help of the parameters (shown in the question marks of the forms), the parameters of
the basic (non-expert) mode and the conditions under which a parameter or category is used at all.

HELP:            'Section.Category.Variable' -> markdown text (what it does, effect of the value, typical values)
BASIC:           parameters shown when the expert mode is off (empty: all parameters are basic)
ACTIVE:          'Section.Category.Variable' -> function(cfg) that is False if the parameter has no effect with the
                 current settings (hidden in the basic mode, marked as inactive in the expert mode)
RELEVANT:        section -> {category: function(cfg)} False if the category is not used by the selected mode
CATEGORY_TITLES: category -> title of its expander
"""


#######################################################################################################################
# Helper for the rules
#######################################################################################################################
def _mode(c):
    return int(c['Exp'].get('type', 0))


def _topo(c):
    return str(c['Top'].get('sourceType', 'B6'))


def _dab(c):
    return _topo(c) == 'DAB'


#######################################################################################################################
# Help texts
#######################################################################################################################
HELP = {
    # ==================================================================================================================
    # Exp: simulation settings
    # ==================================================================================================================
    'Exp.Exp.therFeed':
        "Coupling of the thermal and the electrical model. **1:** the temperatures are fed back into the "
        "temperature-dependent device parameters, so losses and temperatures are solved together. **0:** no "
        "feedback of the temperatures into the electrical model (faster, less accurate at high temperature rise).",
    'Exp.Exp.freqPar':
        "How often the electrical and loss values are updated in the coupled simulation: after every switching cycle "
        "(**fs**, finer, slower) or after every fundamental cycle (**fel**).",
    'Exp.Exp.freqAvg':
        "Optional moving average of the electrical and loss waveforms (src/general/calcAvg.py): **none** keeps the "
        "switching ripple, **fs** averages over one switching period, **fel** over one fundamental period (shows the "
        "mean loss level).",
    'Exp.Exp.fsim':
        "Sampling frequency of the time-domain simulation. It must resolve the switching events: use about "
        "100 × the switching frequency (Par › PWM › Switching Freq.) or more for short dead times. Memory and run "
        "time grow linearly (the sanity check warns about RAM). The templates use 500 kHz for the inverter at about "
        "1 kHz switching and 100 MHz for the DAB at 100 kHz.",
    'Exp.Exp.tol':
        "Convergence tolerance of the steady-state analysis: the iteration stops when the change with respect to the "
        "previous converged result is below this value. Smaller values are more accurate but need more iterations.",
    'Exp.Exp.eps':
        "Numerical zero used to avoid divisions by zero; keep it far below every physical quantity (default 1e-12).",
    'Exp.Exp.int':
        "Number of grid points used to interpolate the device characteristics of the parameter workbooks "
        "(src/data/loadPara*.py). More points give smoother models but slower loading.",

    # ==================================================================================================================
    # Dat: operating point
    # ==================================================================================================================
    'Dat.stat.cyc':
        "Number of fundamental cycles simulated in the stationary analyses; the first cycles contain the start-up "
        "transient. At least 4 (the sanity check stops below 3). More cycles help with large L/R time constants but "
        "increase memory (fsim × cyc / fel samples).",
    'Dat.stat.W':
        "Number of operating points of the sweep analysis. More points give smoother sweep curves; the run time grows "
        "linearly.",
    'Dat.stat.Tj':
        "Junction (core) temperature of all components in the stationary analysis, used for the temperature-"
        "dependent device parameters.",
    'Dat.stat.Tc':
        "Reference (case) temperature of all components in the stationary analysis.",
    'Dat.stat.Po':
        "Set point of the output active power when the controlled quantity is **P** (Setup › Controlled quantity).",
    'Dat.stat.Qo':
        "Set point of the output reactive power when the controlled quantity is **Q**.",
    'Dat.stat.Vo':
        "Set point of the RMS output phase voltage when the controlled quantity is **V**.",
    'Dat.stat.Io':
        "Set point of the RMS output phase current when the controlled quantity is **I**.",
    'Dat.stat.PhiDAB':
        "Phase shift φ between the primary and secondary bridge voltage of the DAB. With single phase shift "
        "modulation the transferred power is P = N·Vdc·Vo·(φ − φ²/π)/(ω·Lk) (src/general/genLoadInput.py), with the "
        "maximum at 90°.",
    'Dat.stat.Mi':
        "Modulation index of the inverter. Above the linear range (Mi > 1 for SPWM, about 1.15 for SVPWM) the output "
        "voltage saturates and low-order harmonics rise. In the sweep analysis the modulation index is swept up to "
        "this value.",
    'Dat.stat.Vdc':
        "DC-link voltage (input voltage of the DAB). Scales the output voltage, the device voltages and the switching "
        "losses.",
    'Dat.stat.phi':
        "Phase angle of the output voltage reference.",
    'Dat.trans.tmax':
        "Simulated time of the transient and closed-loop analyses. Choose several thermal time constants of the "
        "devices to see the temperatures settle; memory grows with fsim × tmax.",
    'Dat.trans.Tc':
        "Reference (case) temperature of all components in the transient analysis.",
    'Dat.trans.Tj':
        "Junction temperature of all components at t = 0 of the transient analysis.",

    # ==================================================================================================================
    # Top: topology and load
    # ==================================================================================================================
    'Top.Top.sourceType':
        "Converter topology: **B2** half bridge, **B4** full bridge (single phase), **B6** two-level three-phase "
        "inverter, **DAB** dual active bridge DC-DC converter (uses the DCDC device files SwiPri/SwiSec/CapPri/"
        "CapSec/Trafo of the Setup page). Load the config *default_DCDC* as starting point for a DAB.",
    'Top.Top.inpFilter': "Activates the RLC filter at the DC input (Rinp, Linp, Cinp).",
    'Top.Top.Rinp': "Resistance of the input filter.",
    'Top.Top.Linp': "Inductance of the input filter.",
    'Top.Top.Cinp': "Capacitance of the input filter.",
    'Top.Top.outFilter': "Activates the RLC filter at the AC output (Rout, Lout, Cout) between converter and load.",
    'Top.Top.Rout': "Resistance of the output filter.",
    'Top.Top.Lout': "Inductance of the output filter.",
    'Top.Top.Cout': "Capacitance of the output filter.",
    'Top.Top.R':
        "Load resistance per phase. Together with L it sets the load angle φ = atan(ωL/R) and the current amplitude.",
    'Top.Top.L':
        "Load inductance per phase. Larger values smooth the current (lower current THD) and increase the load angle.",
    'Top.Top.E': "Amplitude of an induced voltage (back EMF) in the load, e.g. of an electric machine.",
    'Top.Top.phiE': "Phase angle of the induced voltage.",
    'Top.Top.wave':
        "Shape of the reference waveform (src/pwm/genWaveform.py): sinusoidal (inverters), constant (DC-DC), "
        "triangular, rectangular or trapezoidal.",
    'Top.Top.fel':
        "Fundamental frequency of the output (for the DAB: frequency of the transformer voltage). With fsim it sets "
        "the number of samples per cycle.",

    # ==================================================================================================================
    # Par: PWM
    # ==================================================================================================================
    'Par.PWM.type':
        "Modulation method: **FF** fundamental frequency, **CB** carrier based, **SV** space vector, **OPP** optimal "
        "pulse pattern, **OSM** optimal subcycle modulation (inverters); **SPS** single phase shift (DAB).",
    'Par.PWM.upd':
        "Update of the reference: **SE** once per carrier period (single edge), **DE** twice per carrier period "
        "(double edge).",
    'Par.PWM.samp':
        "**NS** natural sampling (reference compared continuously with the carrier) or **RS** regular sampling "
        "(reference sampled and held at the update instants, as in a digital controller).",
    'Par.PWM.tri': "Carrier edge trigger: rising edge, falling edge, symmetrical or asymmetrical modulation.",
    'Par.PWM.int': "Interleaving of the carriers when multiple carriers are used.",
    'Par.PWM.td':
        "Dead time between turning off one switch and turning on the other switch of a half bridge. Causes a "
        "current-direction dependent voltage error and low-order harmonics; must be resolved by fsim.",
    'Par.PWM.tmin': "Minimum on/off time of a pulse.",
    'Par.PWM.loss':
        "**0:** ideal lossless switches, **1:** loss models from the device parameter workbooks (conduction and "
        "switching losses, thermal model).",
    'Par.PWM.swloss':
        "Switching loss model: **0** from the switching energies Eon, Eoff, Erec, **1** from the integration of the "
        "device capacitances Ciss, Coss, Crss.",
    'Par.PWM.sw': "**0:** hard switching. **1:** soft switching (to be implemented).",
    'Par.PWM.fs':
        "Switching frequency. Higher values reduce the current ripple and distortion but increase the switching "
        "losses about linearly; fsim must be about 100 × higher.",
    'Par.PWM.seq':
        "Order of the switching states within a sampling period of the B6 space vector modulation (0 and 7 are the "
        "zero vectors). Changes the number of commutations and the harmonic spectrum.",
    'Par.PWM.zero':
        "Zero-sequence / modulation method of the B6 bridge: SPWM, SVPWM, third harmonic injection (THIPWM4/6) or "
        "discontinuous PWM (DPWM0–3, DPWMMIN/MAX: one phase is clamped per sector, fewer commutations and switching "
        "losses, higher distortion).",

    # ==================================================================================================================
    # Par: control (closed loop)
    # ==================================================================================================================
    'Par.Cont.type':
        "Controller of the closed-loop analysis: **PI** controller or **HY** hysteresis control. In the inverter "
        "topologies (B2, B4, B6) both settings use the hysteresis current control (src/cont/conHys.py, called in "
        "calcCON of the topology classes) — the PI controller is implemented for the DAB only, where it sets the phase "
        "shift.",
    'Par.Cont.fc':
        "Execution frequency of the controller: the measured current is compared with the reference once per control "
        "period 1/fc and the switching state is held in between. The closed loop simulates tmax·fc control periods one "
        "after another, so the run time grows with fc.",
    'Par.Cont.satMin': "Lower limit of the controller output (DAB PI: phase shift in rad, limited to ±π/2).",
    'Par.Cont.satMax': "Upper limit of the controller output (DAB PI: phase shift in rad, limited to ±π/2).",
    'Par.Cont.Kp': "Proportional gain of the PI controller (DAB only): phase shift [rad] per ampere of current error.",
    'Par.Cont.Ki': "Integral gain of the PI controller (DAB only): phase shift [rad] per ampere-second of integrated "
                   "current error; anti-windup stops the integration while the output is saturated.",
    'Par.Cont.hys':
        "Width of the hysteresis band in % of the peak reference current (used as hys/100): the leg switches to −1 "
        "when the current exceeds the reference by more than the band and to +1 when it is lower by more than the "
        "band. Smaller bands give a lower current ripple but a higher, variable switching frequency (limited by fc).",

    # ==================================================================================================================
    # Par: electrical models
    # ==================================================================================================================
    'Par.Elec.SwiMdl':
        "Model of the switch characteristics from the parameter workbook: **con** constant values, **pwl** piecewise "
        "linear, **tab** tabulated (interpolated, most detailed).",
    'Par.Elec.SwiType': "Device technology of the switches (IGBT with antiparallel diode, or MOSFET).",
    'Par.Elec.SwiRecCon':
        "Reverse conduction path: diode only (**D**) or shared by diode and transistor (**DT**, to be implemented).",
    'Par.Elec.SwiRecMdl': "Reverse conduction model: **0** ideal, **1** including the blanking (dead) time.",
    'Par.Elec.SwiPara': "Number of devices in parallel per switch.",
    'Par.Elec.SwiSeries': "Number of devices in series per switch.",
    'Par.Elec.TraMdl': "Model of the transformer parameters (DAB): constant, piecewise linear or tabulated.",
    'Par.Elec.CapMdl': "Model of the capacitor parameters: constant, piecewise linear or tabulated.",
    'Par.Elec.CapType': "Capacitor technology (Elco; MLCC and film to be implemented).",
    'Par.Elec.CapPara': "Number of capacitors in parallel.",
    'Par.Elec.CapSeries': "Number of capacitors in series.",

    # ==================================================================================================================
    # Par: thermal models
    # ==================================================================================================================
    'Par.Ther.Heatsink':
        "**1:** the RC model of the heatsink (thermal resistances and capacities) is used; **0:** the heatsink is "
        "neglected.",
    'Par.Ther.Coupling':
        "Thermal coupling: **0** none, **1** between diode and transistor of a switch, **2** via the whole converter.",
}

#######################################################################################################################
# Basic parameters (expert mode off)
#######################################################################################################################
BASIC = {
    'Exp.Exp.fsim', 'Exp.Exp.therFeed',
    'Dat.stat.cyc', 'Dat.stat.W', 'Dat.stat.Tj', 'Dat.stat.Tc', 'Dat.stat.Po', 'Dat.stat.Qo', 'Dat.stat.Vo',
    'Dat.stat.Io', 'Dat.stat.PhiDAB', 'Dat.stat.Mi', 'Dat.stat.Vdc',
    'Dat.trans.tmax', 'Dat.trans.Tc', 'Dat.trans.Tj',
    'Top.Top.sourceType', 'Top.Top.inpFilter', 'Top.Top.outFilter', 'Top.Top.R', 'Top.Top.L', 'Top.Top.E',
    'Top.Top.wave', 'Top.Top.fel',
    'Par.PWM.type', 'Par.PWM.zero', 'Par.PWM.fs', 'Par.PWM.td', 'Par.PWM.loss',
    'Par.Cont.type', 'Par.Cont.Kp', 'Par.Cont.Ki', 'Par.Cont.hys',
    'Par.Elec.SwiPara', 'Par.Ther.Heatsink', 'Par.Ther.Coupling',
}

#######################################################################################################################
# Activity rules
#######################################################################################################################
ACTIVE = {
    'Exp.Exp.tol': lambda c: _mode(c) == 1,
    'Dat.stat.W': lambda c: _mode(c) == 0,
    'Dat.stat.Po': lambda c: c['Exp'].get('output') == 'P',
    'Dat.stat.Qo': lambda c: c['Exp'].get('output') == 'Q',
    'Dat.stat.Vo': lambda c: c['Exp'].get('output') == 'V',
    'Dat.stat.Io': lambda c: c['Exp'].get('output') == 'I',
    'Dat.stat.PhiDAB': lambda c: _dab(c),
    'Dat.stat.Mi': lambda c: not _dab(c),
    'Top.Top.Rinp': lambda c: c['Top']['inpFilter'] == 1,
    'Top.Top.Linp': lambda c: c['Top']['inpFilter'] == 1,
    'Top.Top.Cinp': lambda c: c['Top']['inpFilter'] == 1,
    'Top.Top.Rout': lambda c: c['Top']['outFilter'] == 1,
    'Top.Top.Lout': lambda c: c['Top']['outFilter'] == 1,
    'Top.Top.Cout': lambda c: c['Top']['outFilter'] == 1,
    'Top.Top.phiE': lambda c: float(c['Top'].get('E') or 0) != 0,
    'Par.PWM.upd': lambda c: not _dab(c),
    'Par.PWM.samp': lambda c: not _dab(c),
    'Par.PWM.tri': lambda c: not _dab(c),
    'Par.PWM.seq': lambda c: _topo(c) == 'B6',
    'Par.PWM.zero': lambda c: _topo(c) == 'B6',
    'Par.PWM.swloss': lambda c: c['Par']['PWM']['loss'] == 1,
    # the PI controller exists for the DAB only, B2/B4/B6 use the hysteresis control for both types (calcCON)
    'Par.Cont.Kp': lambda c: c['Par']['Cont']['type'] == 'PI' and _topo(c) == 'DAB',
    'Par.Cont.Ki': lambda c: c['Par']['Cont']['type'] == 'PI' and _topo(c) == 'DAB',
    'Par.Cont.satMin': lambda c: c['Par']['Cont']['type'] == 'PI' and _topo(c) == 'DAB',
    'Par.Cont.satMax': lambda c: c['Par']['Cont']['type'] == 'PI' and _topo(c) == 'DAB',
    'Par.Cont.hys': lambda c: c['Par']['Cont']['type'] == 'HY' or _topo(c) != 'DAB',
    'Par.Elec.TraMdl': lambda c: _dab(c),
}

#######################################################################################################################
# Categories
#######################################################################################################################
RELEVANT = {
    'Par': {'Cont': lambda c: _mode(c) == 3},
    'Dat': {'trans': lambda c: _mode(c) in (2, 3)},
}
CATEGORY_TITLES = {'PWM': 'PWM', 'Cont': 'Control (closed loop)', 'Elec': 'Electrical models',
                   'Ther': 'Thermal models', 'stat': 'Operating point', 'trans': 'Transient', 'Top': 'Topology & load',
                   'Exp': 'Simulation'}


#######################################################################################################################
# Helper
#######################################################################################################################
def key(section, p):
    return f"{section}.{p['cat']}.{p['var']}"


def is_basic(section, p):
    return not BASIC or key(section, p) in BASIC


def is_active(section, p, cfg):
    rule = ACTIVE.get(key(section, p))
    if rule is None:
        return True
    try:
        return bool(rule(cfg))
    except (KeyError, TypeError, ValueError):
        return True
