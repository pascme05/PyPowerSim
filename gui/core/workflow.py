#######################################################################################################################
#######################################################################################################################
# Title:        Power Electronics Simulation Toolkit (PyPowerSim)
# Topic:        Power Electronics
# File:         workflow.py
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
This function draws the workflow of every operating mode (and sub-mode) as a flow chart (Graphviz, rendered in the
browser) and explains the steps. The nodes show the settings of the current configuration that control the step;
steps that are skipped with the current settings are drawn dashed. Add one function per mode to WORKFLOWS.
"""

#######################################################################################################################
# Import libs
#######################################################################################################################
from gui import adapter as A
from gui.core import forms

import streamlit as st

#######################################################################################################################
# Style
#######################################################################################################################
COL = {'inp': ('#e8f1fc', '#2a78d6'), 'step': ('#ffffff', '#52514e'), 'out': ('#e7f6ef', '#1baf7a'),
       'opt': ('#fdf1e8', '#eb6834'), 'off': ('#f6f6f3', '#c3c2b7')}


def opt(sec, cat, var):
    """'value) label' of a parameter with options."""
    return f"{forms.value(sec, cat, var)}) {forms.option_label(sec, cat, var).split(' (')[0]}"


def val(sec, cat, var, unit=''):
    v = forms.value(sec, cat, var)
    if isinstance(v, float) and forms.find(sec, cat, var)['percent']:
        return f"{100 * v:g} %"
    return f"{v:g}{(' ' + unit) if unit else ''}" if isinstance(v, (int, float)) else str(v)


class Flow:
    """Small builder of a Graphviz flow chart."""

    def __init__(self, rankdir='TB'):
        self.nodes, self.edges, self.ranks = [], [], []
        self.rankdir = rankdir

    def node(self, nid, title, lines=(), kind='step', active=True):
        fill, line = COL[kind if active else 'off']
        txt = f"<<B>{_esc(title)}</B>" + ''.join(f"<BR/><FONT POINT-SIZE='10' COLOR='#52514e'>{_esc(l)}</FONT>"
                                                  for l in lines if l) + ">"
        shape = {'inp': 'note', 'out': 'box', 'opt': 'box', 'step': 'box'}[kind]
        style = 'filled,rounded' + (',dashed' if not active else '') + (',bold' if kind == 'out' else '')
        self.nodes.append(f'"{nid}" [label={txt}, shape={shape}, style="{style}", fillcolor="{fill}", '
                          f'color="{line}", fontcolor="{"#0b0b0b" if active else "#898781"}"];')
        return self

    def edge(self, a, b, label='', active=True):
        lab = f', label="{label}", fontsize=9, fontcolor="#52514e"' if label else ''
        self.edges.append(f'"{a}" -> "{b}" [color="{"#52514e" if active else "#c3c2b7"}"'
                          f'{", style=dashed" if not active else ""}{lab}];')
        return self

    def same(self, *ids):
        self.ranks.append('{rank=same; ' + ' '.join(f'"{i}";' for i in ids) + '}')
        return self

    def dot(self):
        return ('digraph G {\n'
                f'  rankdir={self.rankdir}; bgcolor="transparent"; nodesep=0.35; ranksep=0.4;\n'
                '  node [fontname="Helvetica", fontsize=12, margin="0.18,0.08"];\n'
                '  edge [arrowsize=0.7, penwidth=1.2];\n  ' + '\n  '.join(self.nodes + self.edges + self.ranks) +
                '\n}')


def _esc(s):
    return str(s).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


def _legend():
    st.caption(":blue[■] inputs (workbooks, config) · ■ processing step · :orange[■] option / loop · "
               ":green[■] results · dashed: skipped with the current settings")


#######################################################################################################################
# Workflows: mode -> function(cfg, sub) returning (Flow, [(step title, explanation), ...])
#######################################################################################################################
def _example(cfg, sub):
    f = Flow()
    f.node('cfg', 'Configuration', [A.NAME], kind='inp')
    f.node('sim', 'Simulation', ['main(setup, path)'])
    f.node('res', 'Results', ['figures, data'], kind='out')
    f.edge('cfg', 'sim').edge('sim', 'res')
    return f, [('Configuration', "parameters of the GUI"), ('Simulation', "toolkit"), ('Results', "outputs")]


#######################################################################################################################
# PyPowerSim (src/main.py)
#######################################################################################################################
CALC = {0: ('Sweep', 'calcSweep', "the operating point is swept over W points (modulation index or phase shift); "
                                  "for every point the stationary waveforms and the analytical solution are "
                                  "calculated"),
        1: ('Steady state', 'calcSteady', "cyc fundamental periods are simulated and iterated until the result "
                                          "changes less than the tolerance; the last period is the stationary result"),
        2: ('Transient', 'calcTrans', "the waveforms are simulated over tmax with the thermal RC networks starting "
                                      "at Tj, so the heating of the devices is resolved"),
        3: ('Closed loop', 'calcClose', "the output is controlled by the PI or hysteresis controller (Par › Cont) "
                                        "over tmax")}


def _pps(cfg, sub):
    m = int(cfg['Exp']['type'])
    E, T, P = cfg['Exp'], cfg['Top'], cfg['Par']
    dab = T['sourceType'] == 'DAB'
    loss = P['PWM']['loss'] == 1
    title, fnc, txt = CALC[m]
    files = [E['SwiPri'], E['SwiSec'], E['CapPri'], E['CapSec'], E['Trafo']] if dab else [E['Swi'], E['Cap']]

    f = Flow()
    f.node('cfg', 'Configuration', [f"GUI → config/_gui_run.xlsx", f"fsim = {val('Exp', 'Exp', 'fsim', 'Hz')}"],
           kind='inp')
    f.node('para', 'Device parameters', [', '.join(files)], kind='inp')
    f.node('load', 'Loading', ['loadSetup, loadPara', f"model: {opt('Par', 'Elec', 'SwiMdl')}"])
    f.node('topo', 'Topology', [opt('Top', 'Top', 'sourceType'), f"load R = {val('Top', 'Top', 'R', 'Ω')}, "
                                                                 f"L = {val('Top', 'Top', 'L', 'H')}"])
    f.node('op', 'Operating point', [f"controlled: {E['output']}",
                                     f"PhiDAB = {val('Dat', 'stat', 'PhiDAB', '°')}" if dab else
                                     f"Mi = {val('Dat', 'stat', 'Mi')}", f"Vdc = {val('Dat', 'stat', 'Vdc', 'V')}"])
    f.node('pwm', 'PWM', [opt('Par', 'PWM', 'type'), f"fs = {val('Par', 'PWM', 'fs', 'Hz')}"] +
           ([] if dab or T['sourceType'] != 'B6' else [P['PWM']['zero']]))
    f.node('calc', title, [f"{fnc}{'_DCDC' if dab else ''}"] +
           ([f"W = {val('Dat', 'stat', 'W')} points"] if m == 0 else
            [f"cyc = {val('Dat', 'stat', 'cyc')}"] if m == 1 else [f"tmax = {val('Dat', 'trans', 'tmax', 's')}"]) +
           ([opt('Par', 'Cont', 'type')] if m == 3 else []), kind='opt')
    f.node('el', 'Electrical & losses', [opt('Par', 'PWM', 'loss'), opt('Par', 'PWM', 'swloss') if loss else ''])
    f.node('th', 'Thermal', [opt('Par', 'Ther', 'Heatsink').split(' (')[0][:40], opt('Par', 'Ther', 'Coupling')[:40]],
           active=loss)
    f.node('out', 'Results', ['waveforms, spectra, losses, temperatures'] + (['sweep curves'] if m == 0 else []),
           kind='out')
    f.node('post', 'Outputs', [f"plots: {E['plot']}", f"save: {E['save']}"], kind='out', active=bool(E['plot'] or
                                                                                                     E['save']))
    f.edge('cfg', 'load').edge('para', 'load').same('cfg', 'para')
    f.edge('load', 'topo').edge('topo', 'op').edge('op', 'pwm').edge('pwm', 'calc').edge('calc', 'el')
    f.edge('el', 'th', active=loss)
    f.edge('th', 'el', label='therFeed' if E['therFeed'] == 1 else '', active=loss and E['therFeed'] == 1)
    f.edge('el', 'out').edge('out', 'post', active=bool(E['plot'] or E['save']))
    steps = [('Configuration', "the parameters of the GUI are written into a copy of the config workbook, which "
                               "the toolkit loads like a start script"),
             ('Device parameters', "switches, capacitors" + (" and transformer" if dab else "") +
              " from para/ (selected on the Setup page)"),
             ('Topology', "the converter class (B2/B4/B6/DAB) is initialised with the load and filters"),
             ('Operating point', "the set point of the controlled quantity is converted into the modulation index or "
                                 "phase shift (src/general/genLoadInput.py)"),
             ('PWM', "switching functions from the reference and the carrier"),
             (title, txt),
             ('Electrical & losses', "device currents/voltages, conduction and switching losses" if loss else
              "ideal switches (Losses = 0): no losses and no thermal model"),
             ('Thermal', "Foster RC networks of the devices (and heatsink); with thermal coupling the temperatures "
                         "are fed back into the device parameters"),
             ('Results', "shown interactively in Results; figures of the toolkit under Figures & files")]
    return f, steps


WORKFLOWS = {m: _pps for m in CALC}


#######################################################################################################################
# View
#######################################################################################################################
def view(m):
    cfg = forms.cfg()
    sub = None
    spec = A.MODES[m].get('sub')
    if spec:
        p = forms.find(*spec)
        cur = forms.value(*spec)
        if p['options']:
            names = p['options']
            sub = st.segmented_control("Sub-mode", list(names), default=cur, key=f"wf|{m}",
                                       format_func=lambda k: f"{names[k]}" + (' (selected)' if k == cur else ''))
            sub = cur if sub is None else sub
            if sub != cur:
                st.caption(f"Showing sub-mode {sub}; the configuration runs sub-mode {cur} — change it in "
                           f"*Parameters*.")
    f, steps = WORKFLOWS[m](cfg, sub)
    c1, c2 = st.columns([3, 2])
    with c1:
        st.graphviz_chart(f.dot(), width='content')
        _legend()
    with c2:
        st.markdown("##### What happens")
        for i, (t, txt) in enumerate(steps, start=1):
            st.markdown(f"**{i}. {t}** — {txt}")
        st.caption("The boxes show the settings of the current configuration; change them in *Parameters*, the "
                   "diagram follows.")
