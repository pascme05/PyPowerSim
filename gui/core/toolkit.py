#######################################################################################################################
#######################################################################################################################
# Title:        Power Electronics Simulation Toolkit (PyPowerSim)
# Topic:        Power Electronics
# File:         toolkit.py
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
Access to the functions of the toolkit for the theory explorers and input charts: the pre-processing chain of
src/main.py (loadPara, sanityInput, genTF, genLoadInput, initTopo) started from the GUI configuration, and the PWM and
time-domain functions of the topology classes for one fundamental period. The console output of the toolkit is
suppressed; results are cached per configuration.
"""

#######################################################################################################################
# Import libs
#######################################################################################################################
from gui import adapter as A

import contextlib
import copy
import io
import json
import math

import numpy as np
import streamlit as st


#######################################################################################################################
# Pre-processing
#######################################################################################################################
def _quiet():
    return contextlib.redirect_stdout(io.StringIO())


@st.cache_resource(show_spinner=False, max_entries=8)
def _prepare(cfg_json):
    """(setup, para, mdl, top) after the pre-processing of src/main.py for one fundamental period (cyc = 1)."""
    from src.data.loadPara import loadPara
    from src.general.sanityCheck import sanityInput
    from src.general.genTF import genTF
    from src.general.genLoadInput import genLoadInput
    from src.topo.initTopo import initTopo
    s = json.loads(cfg_json)
    s['Dat']['stat']['cyc'] = 1
    with _quiet():
        para = loadPara(s, A.paths())
        s = sanityInput(para, s)
        mdl = genTF(para, s)
        s = genLoadInput(s, para)
        top = initTopo(s['Top']['sourceType'], s, para)
    return s, para, mdl, top


def prepare(cfg, **changes):
    """
    Pre-processing for the configuration cfg with changes {'Sec.cat.var': value} (e.g. {'Exp.fsim': 2e5}). The
    returned objects are shared by the cache — copy before modifying them.
    """
    s = copy.deepcopy(cfg)
    for k, v in changes.items():
        d = s
        *path, last = k.split('.')
        for p in path:
            d = d[p]
        d[last] = v
    return _prepare(json.dumps(s, sort_keys=True, default=str))


#######################################################################################################################
# One fundamental period with the toolkit
#######################################################################################################################
@st.cache_data(show_spinner=False, max_entries=32)
def _period(cfg_json, Mi):
    s, _, mdl, top = _prepare(cfg_json)
    s = copy.deepcopy(s)
    fel, fsim = s['Top']['fel'], s['Exp']['fsim']
    N = int(fsim / fel)
    t = np.linspace(0, 1 / fel, N + 1)
    with _quiet():
        v_ref, e_ref, _ = top.calcRef(s['Top']['E'], math.radians(s['Top']['phiE']),
                                      math.radians(s['Dat']['stat']['phi']), [], s)
        xs, xsh, sw, c, x, xN0 = top.calcPWM(v_ref, t, Mi, s)
        ac, dc, _ = top.calcTime(sw, e_ref, t, Mi, mdl, 0, len(t), [], 0, s)
    keys = [k for k in ('A', 'B', 'C') if k in sw]
    return {'t': t, 'c': np.asarray(c, dtype=float), 'x': {k: np.asarray(x[k], dtype=float) for k in keys},
            's': {k: np.asarray(sw[k], dtype=float) for k in keys}, 'n0': np.asarray(xN0, dtype=float),
            'ac': {k: np.asarray(v, dtype=float) for k, v in ac.items() if not isinstance(v, dict)},
            'dc': {k: np.asarray(v, dtype=float) for k, v in dc.items() if not isinstance(v, dict)},
            'Vdc': float(np.abs(s['Dat']['stat']['Vdc'])), 'fel': fel, 'fsim': fsim}


def period(cfg, Mi, **changes):
    """
    Reference, carrier, switching functions and ideal output (calcRef, calcPWM, calcTime of the topology class) of
    one fundamental period for the modulation index Mi; inverter topologies only.
    """
    s = copy.deepcopy(cfg)
    for k, v in changes.items():
        d = s
        *path, last = k.split('.')
        for p in path:
            d = d[p]
        d[last] = v
    return _period(json.dumps(s, sort_keys=True, default=str), float(Mi))
