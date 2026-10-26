#######################################################################################################################
#######################################################################################################################
# Title:        Power Electronics Simulation Toolkit (PyPowerSim)
# Topic:        Power Electronics
# File:         adapter.py
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
Connection of the GUI to the toolkit. This is the only file (together with gui/core/help.py) that is specific to the
toolkit: names, config workbook sheets, run settings of start_INV.py/start_DCDC.py, operating modes, parameter folders
and the call of the simulation. All other GUI modules are generic and read their settings from here.
"""

#######################################################################################################################
# Import libs
#######################################################################################################################
import contextlib
import copy
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#######################################################################################################################
# Names
#######################################################################################################################
NAME = 'PyPowerSim'                                                   # short name (page title, JSON meta)
TITLE = 'Power Electronics Simulation Toolkit'                        # full name (Setup page)
TAGLINE = "Simulation of inverters (B2, B4, B6) and DC-DC converters (DAB): PWM, waveforms, losses, temperatures and " \
          "spectra. Set up the run here, then open an operating mode in the sidebar. All parameters are edited in the " \
          "GUI and stored as JSON config (sidebar); the device parameters are read from the workbooks in para/."
SRC_PATHS = ['', 'src']                                               # added to sys.path (relative to ROOT)
LOGO = None                                                           # e.g. 'docu/logo.png' (optional)
DOC = 'docu/GUI.md'                                                   # user guide shown on the Documentation page

#######################################################################################################################
# Config workbook (Name | Category | Description | Variable | Value | Unit)
#######################################################################################################################
CONFIG_DIR = 'config'                                                 # folder of the config workbooks and JSON
CONFIG_WORKBOOK = 'config/default.xlsx'                               # template, source of schema and defaults
# sheet -> (section of setup, nested): same structure as src/data/loadSetup.py
SHEETS = {
    'Exp': ('Exp', False),                                            # setup['Exp'][var]
    'Dat': ('Dat', True),                                             # setup['Dat']['stat'|'trans'][var]
    'Top': ('Top', False),                                            # setup['Top'][var]
    'Par': ('Par', True),                                             # setup['Par']['PWM'|'Cont'|'Elec'|'Ther'][var]
}
INTS = {('Exp', 'Exp', 'int'), ('Dat', 'stat', 'cyc'), ('Dat', 'stat', 'W'), ('Par', 'Elec', 'SwiPara'),
        ('Par', 'Elec', 'SwiSeries'), ('Par', 'Elec', 'CapPara'), ('Par', 'Elec', 'CapSeries')}
STRS = set()                                                          # {(section, cat, var)} forced text
NO_OPTIONS = set()                                                    # descriptions enumerate cases, value continuous
OPTIONAL = set()                                                      # empty value allowed
# options used by the code (src/) that the workbook does not list
OPTIONS = {
    ('Top', 'Top', 'sourceType'): {'DAB': 'dual active bridge (DCDC)'},
    ('Par', 'PWM', 'zero'): {'SPWM': 'sinusoidal PWM', 'SVPWM': 'space vector PWM',
                             'THIPWM4': 'third harmonic injection (1/4)', 'THIPWM6': 'third harmonic injection (1/6)',
                             'DPWM0': 'discontinuous PWM 0', 'DPWM1': 'discontinuous PWM 1',
                             'DPWM2': 'discontinuous PWM 2', 'DPWM3': 'discontinuous PWM 3',
                             'DPWMMAX': 'discontinuous PWM (max)', 'DPWMMIN': 'discontinuous PWM (min)'},
    ('Par', 'PWM', 'seq'): {'0127': 'states 0-1-2-7', '0121': 'states 0-1-2-1', '7212': 'states 7-2-1-2',
                            '1012': 'states 1-0-1-2', '2721': 'states 2-7-2-1'},
    ('Par', 'Elec', 'SwiType'): {'IGBT': 'IGBT with antiparallel diode', 'MOSFET': 'MOSFET'},
    # src/pwm/genWaveform.py: sin, con, tri, rec, tra (the workbook names trapezoidal 'tar', which is not handled)
    ('Top', 'Top', 'wave'): {'tar': None, 'rec': 'rectangular', 'tra': 'trapezoidal'},
}
LABELS = {('Par', 'PWM', 'sw'): 'Switching Mode'}
PERCENT_AS_FRACTION = False                                           # hys [%] is used as hys/100 in src/
# further ready-to-load configs from workbooks of the same layout: name -> (workbook, run settings of start_DCDC.py)
EXTRA_CONFIGS = {'default_DCDC': ('config/default_DCDC.xlsx', {'Name': 'DAB_SPS', 'output': 'V', 'type': 1,
                                                               'conf': 'default_DCDC', 'plot': 3})}

#######################################################################################################################
# Run settings (formerly start_INV.py / start_DCDC.py)
#######################################################################################################################
RUN_SECTION = 'Exp'                                                   # section of setup holding the run settings
CONFIG_VAR = 'conf'                                                   # run setting naming the config workbook
MODE_VAR = 'type'                                                     # run setting of the operating mode
RUN = {'Name': 'default', 'Author': 'Pascal Schirmer', 'debug': 0, 'output': 'Mi', 'type': 0, 'conf': 'default',
       'Swi': 'IKQ75N120CS6', 'Cap': 'Elco', 'SwiPri': 'A2F12M12W2', 'SwiSec': 'A2F06M12W2', 'CapPri': 'Cin',
       'CapSec': 'Cout', 'Trafo': 'trafoDAB', 'plot': 2, 'plotGen': 0, 'save': 0}
# (var, kind, label, help): kind 'text', 'toggle', 'int', 'file:<folder key>', 'select:a|b|c', 'expert:<kind>'
RUN_FIELDS = [
    ('Name', 'text', "Simulation name", "Name of the results folder results/<Name>/ (used when saving)."),
    ('Author', 'text', "Author", "Responsible person, written into the saved setup."),
    ('output', 'select:Mi|Phi|V|I|P|Q', "Controlled quantity",
     "Mi: modulation index, Phi: phase shift (DAB), V: output voltage, I: output current, P: active power, "
     "Q: reactive power. The set points are the operating point parameters (Dat › stat: Mi, PhiDAB, Vo, Io, Po, "
     "Qo)."),
    ('plot', 'select:0|1|2|3', "Plots",
     "0: no plots, 1: results, 2: results with analytic solutions, 3: topology specific plots. Figures are shown "
     "under Results › Figures & files."),
    ('plotGen', 'select:0|1', "Model plots", "1: loss and thermal models of the devices are plotted in addition."),
    ('save', 'select:0|1', "Save results", "1: results and setup are written to results/<Name>/."),
    ('debug', 'expert:select:0|1|2', "Debug level", "Amount of console output (0 normal, 2 most detailed)."),
    ('Swi', 'file:swi', "Switches", "Parameter workbook of the switching devices in para/Swi/ (B2, B4, B6)."),
    ('Cap', 'file:cap', "DC-link capacitor", "Parameter workbook of the DC-link capacitor in para/Cap/."),
    ('SwiPri', 'file:swi', "Switches primary bridge", "Switching devices of the primary bridge (DAB) in para/Swi/."),
    ('SwiSec', 'file:swi', "Switches secondary bridge", "Switching devices of the secondary bridge (DAB)."),
    ('CapPri', 'file:cap', "Input capacitor", "Input capacitor of the DAB in para/Cap/."),
    ('CapSec', 'file:cap', "Output capacitor", "Output capacitor of the DAB in para/Cap/."),
    ('Trafo', 'file:tra', "Transformer", "Transformer of the DAB in para/Tra/."),
]
# input files grouped by converter type: title, caption, icon, run settings, active(cfg) (the topology selects the group)
FILE_GROUPS = [
    dict(title='Inverter', text='DC/AC converter: half bridge B2, full bridge B4, three-phase bridge B6',
         icon=':material/electric_meter:', fields=['Swi', 'Cap'],
         active=lambda cfg: cfg['Top']['sourceType'] != 'DAB'),
    dict(title='DC-DC converter', text='Dual active bridge (DAB): two full bridges coupled by a transformer',
         icon=':material/compare_arrows:', fields=['SwiPri', 'SwiSec', 'CapPri', 'CapSec', 'Trafo'],
         active=lambda cfg: cfg['Top']['sourceType'] == 'DAB'),
]
GROUP_SELECTOR = ('Top', 'Top', 'sourceType')                       # parameter that selects the active group

#######################################################################################################################
# Operating modes
#######################################################################################################################
_EXTRA = ['Top', 'Dat']
MODES = {
    0: dict(title='Sweep', section='Par', extra=_EXTRA, icon=':material/stacked_line_chart:',
            text="Sweep of the modulation index (inverter) or the phase shift (DAB): distortion, losses and "
                 "efficiency over the operating range (W points, Dat › stat).", sub=None),
    1: dict(title='Steady state', section='Par', extra=_EXTRA, icon=':material/show_chart:',
            text="Stationary waveforms of one operating point (cyc fundamental cycles): switching functions, "
                 "currents, voltages, losses, temperatures and spectra.", sub=None),
    2: dict(title='Transient', section='Par', extra=_EXTRA, icon=':material/timeline:',
            text="Transient analysis of an operating point incl. the thermal response of the devices "
                 "(Dat › trans: initial temperatures).", sub=None),
    3: dict(title='Closed loop', section='Par', extra=_EXTRA, icon=':material/autorenew:',
            text="Closed-loop control of the output (PI or hysteresis controller, Par › Cont) incl. waveforms "
                 "and losses.", sub=None),
}

#######################################################################################################################
# Files
#######################################################################################################################
FOLDERS = {'swi': 'para/Swi', 'cap': 'para/Cap', 'tra': 'para/Tra'}  # parameter workbooks: key -> folder
FILE_EXT = ('.xlsx',)
RESULTS_DIR = 'results'
RESULT_PROJECTS = False                                               # saveResults writes results/<Name>/


def results_root(cfg):
    """Folder in which a run creates its results folder (new or updated sub-folders are detected after the run)."""
    return os.path.join(ROOT, RESULTS_DIR)


#######################################################################################################################
# Toolkit
#######################################################################################################################
# The loader (src/data/loadSetup.py) is called inside main(): the GUI writes the config into config/_gui_run.xlsx
FEED = 'xlsx'

# Calculation functions called by main(); their outputs are recorded for the interactive results
CALCS = ('calcSweep', 'calcSteady', 'calcTrans', 'calcClose',
         'calcSweep_DCDC', 'calcSteady_DCDC', 'calcTrans_DCDC', 'calcClose_DCDC')


def init_setup():
    """Empty setup of the toolkit."""
    from src.general.helpFnc import initSetup
    return initSetup()


def paths():
    """Path variable of the toolkit; called with ROOT as working directory (initPath uses os.getcwd())."""
    from src.general.helpFnc import initPath
    return initPath(os.path.basename(ROOT))


def cast(setup):
    """Type casting for FEED='dict' (not used, loadSetup reads the workbook)."""
    return setup


def run(setup, path, mode):
    """
    Runs main() of the toolkit. main() returns nothing, so the calculation functions it calls are wrapped for the
    duration of the run and their outputs (time domain, frequency domain, sweep) are returned for the GUI.
    """
    import src.main as M
    rec = {}

    def wrap(name, fnc):
        def inner(*args, **kw):                                       # (top, mdl, para, setup)
            out = fnc(*args, **kw)
            # copies: the plot functions cut the time series to the last period in place
            rec.update({'calc': name, 'setup': copy.deepcopy(args[-1]), 'time': copy.deepcopy(out[0]),
                        'freq': copy.deepcopy(out[1]), 'sweep': copy.deepcopy(out[2]) if len(out) > 2 else None})
            return out
        return inner

    orig = {n: getattr(M, n) for n in CALCS}
    try:
        for n, f in orig.items():
            setattr(M, n, wrap(n, f))
        M.main(setup, path)
    finally:
        for n, f in orig.items():
            setattr(M, n, f)
    return rec or None


MAX_SAMPLES = 40000                                                   # per signal in gui_run.pkl


def compact(data, mode):
    """
    Part of the data stored with a run (gui_run.pkl) to re-open it in the results browser: long signals are thinned
    to at most MAX_SAMPLES samples (the session keeps the full resolution).
    """
    import numpy as np
    import pandas as pd

    def thin(x):
        if isinstance(x, dict):
            return {k: thin(v) for k, v in x.items()}
        if isinstance(x, (np.ndarray, pd.Series, pd.DataFrame)) and len(x) > MAX_SAMPLES:
            s = int(np.ceil(len(x) / MAX_SAMPLES))
            return x[::s] if isinstance(x, np.ndarray) else x.iloc[::s].reset_index(drop=True)
        return x
    return thin(data)


@contextlib.contextmanager
def session():
    """Context of a run (warnings off)."""
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        yield
