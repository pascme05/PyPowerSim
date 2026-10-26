#######################################################################################################################
#######################################################################################################################
# Title:        Power Electronics Simulation Toolkit (PyPowerSim)
# Topic:        Power Electronics
# File:         forms.py
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
This function generates the parameter forms of the GUI from the schema. Every parameter is edited in place in the
session configuration (st.session_state.cfg). Values with unit '%' are fractions in the configuration and shown as
percent in the form. With the expert mode off only the basic parameters that have an effect for the current
settings are shown; hidden parameters keep their value.
"""

#######################################################################################################################
# Import libs
#######################################################################################################################
# ==============================================================================
# Internal
# ==============================================================================
from gui import adapter as A
from gui.core.config import default_config, load_schema, merge, slot
from gui.core.help import HELP, CATEGORY_TITLES, RELEVANT, key as help_key, is_basic, is_active

# ==============================================================================
# External
# ==============================================================================
import math

import streamlit as st

WIDGET_PREFIXES = ('p|', 'g|', 'sb|name', 'sb|pick')


#######################################################################################################################
# Session
#######################################################################################################################
def schema():
    if 'schema' not in st.session_state:
        st.session_state.schema = load_schema()
    return st.session_state.schema


def set_config(cfg, name):
    """Replaces the session configuration; the widgets are re-initialised from it."""
    st.session_state.cfg = merge(default_config(), cfg)
    st.session_state.cfg_name = name
    for k in [k for k in st.session_state if isinstance(k, str) and k.startswith(WIDGET_PREFIXES)]:
        del st.session_state[k]


def cfg():
    if 'cfg' not in st.session_state:
        set_config(default_config(), 'default')
    return st.session_state.cfg


def expert():
    """Expert mode: all parameters are shown (toggle in the sidebar)."""
    return bool(st.session_state.get('expert', False))


def find(section, cat, var):
    for p in schema()[section]:
        if p['cat'] == cat and p['var'] == var:
            return p
    raise KeyError(f"{section}.{cat}.{var}")


#######################################################################################################################
# Widgets
#######################################################################################################################
def _step(v):
    """Step size of a number input in the order of magnitude of the value."""
    if not v:
        return 0.01
    return 10.0 ** math.floor(math.log10(abs(v))) / 10


def _help(section, p, active=True):
    """Question mark of a parameter: extended help, description of the template, variable and default."""
    k = help_key(section, p)
    txt = HELP.get(k) or p['desc'] or p['name']
    if HELP.get(k) and p['desc'] and len(p['desc']) > 60 and p['desc'][:40] not in HELP[k]:
        txt += f"\n\n*Template:* {p['desc']}"
    if p['optional']:
        txt += "\n\nOptional: leave empty if not used."
    if not active:
        txt += "\n\n**Inactive:** this parameter has no effect with the current settings."
    dflt = p['default']
    if p['percent'] and isinstance(dflt, (int, float)):
        dflt = f"{100 * dflt:g} %"
    elif p['options'] and dflt in p['options']:
        dflt = f"{dflt}) {p['options'][dflt]}"
    txt += f"\n\n`{k}` · default: {dflt if dflt not in (None, '') else '(empty)'}"
    return txt.replace('\n', '  \n')


def param(section, p, c=None, disabled=False):
    """Widget of one parameter, the value is written back into the configuration."""
    c = c or cfg()
    store = slot(c, section, p)
    cur = store.get(p['var'], p['default'])
    key = f"p|{section}|{p['cat']}|{p['var']}"
    unit = '%' if p['percent'] else p['unit']
    label = str(p['name']) + (f" [{unit}]" if unit else '')
    active = is_active(section, p, c)
    if not active:
        label += " · inactive"
    hlp = _help(section, p, active)

    # Options
    if p['options']:
        keys = list(p['options'])
        if cur not in keys:
            keys = [cur] + keys
        val = st.selectbox(label, keys, index=keys.index(cur), key=key, help=hlp, disabled=disabled,
                           format_func=lambda k: f"{k}) {p['options'].get(k, '(not defined)')}")

    # Boolean
    elif p['type'] == 'bool':
        val = st.toggle(label, value=bool(cur), key=key, help=hlp, disabled=disabled)

    # Integer
    elif p['type'] == 'int':
        v = None if cur is None else int(round(float(cur)))
        val = st.number_input(label, value=v, step=1, key=key, help=hlp, disabled=disabled)

    # Float
    elif p['type'] == 'float':
        try:
            v = None if cur is None else float(cur) * (100 if p['percent'] else 1)
        except (TypeError, ValueError):
            v = None
        if v is None and not p['optional']:
            v = 0.0
        val = st.number_input(label, value=v, step=_step(v), format='%g', key=key, help=hlp, disabled=disabled)
        if val is not None and p['percent']:
            val = round(val / 100, 12)

    # Text
    else:
        val = st.text_input(label, value='' if cur is None else str(cur), key=key, help=hlp, disabled=disabled)
        val = val.strip() or (None if p['optional'] else '')

    store[p['var']] = val
    return val


def visible(name, p, c=None):
    """Parameter shown in the form (all in the expert mode, else the basic ones that have an effect)."""
    c = c or cfg()
    return expert() or (is_basic(name, p) and is_active(name, p, c))


def section(name, cats=None, cols=3, exclude=()):
    """
    Form of a config section grouped by category. Categories that are not used by the selected sub-mode are collapsed
    and marked (expert mode) or hidden (basic mode).
    """
    c = cfg()
    if name not in schema():
        st.info(f"No parameters for section '{name}' in the schema (adapter.SHEETS).")
        return
    params = [p for p in schema()[name] if (cats is None or p['cat'] in cats) and p['var'] not in exclude
              and not (name == A.RUN_SECTION and p['var'] in A.RUN)]
    shown = [p for p in params if visible(name, p, c)]
    hidden = len(params) - len(shown)
    for cat in dict.fromkeys(p['cat'] for p in params):
        rel = RELEVANT.get(name, {}).get(cat)
        used = True if rel is None else bool(rel(c))
        title = CATEGORY_TITLES.get(cat, cat)
        cat_params = [p for p in shown if p['cat'] == cat]
        if not cat_params or (not used and not expert()):
            hidden += len(cat_params)
            continue
        with st.expander(f"**{title}**" + ("" if used else "  ·  _not used by the selected sub-mode_"),
                         expanded=used):
            grid = st.columns(cols)
            for i, p in enumerate(cat_params):
                with grid[i % cols]:
                    param(name, p, c)
    if hidden and not expert():
        st.caption(f":material/visibility_off: {hidden} expert or inactive parameters hidden — switch on "
                   f"**Expert mode** in the sidebar to edit them (their current values are still used).")


def reset_section(name, label='Reset to defaults'):
    """Button restoring the defaults of a section (run settings are kept)."""
    if st.button(label, key=f"reset|{name}", icon=':material/restart_alt:', type='tertiary'):
        d = default_config().get(name, {})
        if name == A.RUN_SECTION:
            d = {**d, **{k: v for k, v in st.session_state.cfg[name].items() if k in A.RUN}}
        st.session_state.cfg[name] = d
        for k in [k for k in st.session_state if isinstance(k, str) and k.startswith(f"p|{name}|")]:
            del st.session_state[k]
        st.rerun()


def value(section, cat, var):
    return slot(cfg(), section, find(section, cat, var))[var]


def option_label(section, cat, var):
    """Label of the selected option of a parameter, e.g. for summaries and workflow charts."""
    p = find(section, cat, var)
    v = value(section, cat, var)
    return p['options'].get(v, str(v)) if p['options'] else str(v)
