#######################################################################################################################
#######################################################################################################################
# Title:        Power Electronics Simulation Toolkit (PyPowerSim)
# Topic:        Power Electronics
# File:         viz_results.py
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
This function shows the results of a run: interactive charts (RENDER: mode -> function(data, setup, key), generic key
figures otherwise), the figures of the toolkit (captured from matplotlib and found in the results folder), the files
of the results folder and the console log.
"""

#######################################################################################################################
# Import libs
#######################################################################################################################
from gui.core import files
from gui.core.style import fmt, table, show, line, color, downsample, stem

import os

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st


#######################################################################################################################
# Helper
#######################################################################################################################
def arr(x):
    return np.asarray(pd.to_numeric(pd.Series(np.ravel(np.asarray(x, dtype=object))), errors='coerce'), dtype=float)


def scalars(d, prefix='', depth=0):
    """Flat table of the scalar entries of a nested result dictionary."""
    rows = []
    if depth > 6:
        return rows
    for k, v in (d or {}).items() if isinstance(d, dict) else []:
        if isinstance(v, dict):
            rows += scalars(v, f"{prefix}{k}.", depth + 1)
        elif isinstance(v, (int, float, np.integer, np.floating, str)) and not isinstance(v, (bool, np.bool_)):
            rows.append({'quantity': f"{prefix}{k}", 'value': v if isinstance(v, str) else float(v)})
    return rows


def frames(d, prefix='', depth=0):
    """DataFrames and Series of a nested result dictionary {path: DataFrame}."""
    out = {}
    if depth > 6 or not isinstance(d, dict):
        return out
    for k, v in d.items():
        if isinstance(v, pd.DataFrame) and not v.empty:
            out[f"{prefix}{k}"] = v
        elif isinstance(v, pd.Series) and not v.empty:
            out[f"{prefix}{k}"] = v.to_frame()
        elif isinstance(v, dict):
            out.update(frames(v, f"{prefix}{k}.", depth + 1))
    return out


def generic(data, setup, key):
    """Key figures and tables of a result dictionary (until a dedicated view is added to RENDER)."""
    rows = scalars(data)
    if rows:
        st.markdown("**Key figures**")
        table(pd.DataFrame(rows), hide_index=True, height=min(420, 38 + 35 * len(rows)))
    tabs = frames(data)
    if tabs:
        pick = st.selectbox("Table", list(tabs), key=f"{key}|tab")
        table(tabs[pick], height=360)
    if not rows and not tabs:
        st.info("The run returned no scalar results or tables — see *Figures & files*.")


#######################################################################################################################
# Figures and files
#######################################################################################################################
def gallery(run, key='gal'):
    figs = [(t, b) for t, b in run.get('figs') or []]
    disk, reports, datas = [], [], []
    for f in run.get('folders') or []:
        rf = files.run_files(f)
        disk += rf['figures']
        reports += rf['report']
        datas += rf['data']
    items = figs + [(os.path.splitext(os.path.basename(f))[0].replace('_', ' '), f) for f in disk]
    if not items and not reports and not datas:
        st.info("No figures or files for this run (enable rendering/saving on the Setup page).")
        return

    if items:
        st.markdown(f"**Figures of the toolkit** ({len(items)})")
        mode = st.segmented_control("View", ['Grid', 'Single'], default='Grid', key=f"{key}|view",
                                    label_visibility='collapsed')
        if mode == 'Single':
            i = st.selectbox("Figure", range(len(items)), format_func=lambda i: items[i][0], key=f"{key}|one")
            st.image(items[i][1], caption=items[i][0], width='stretch')
        else:
            cols = st.columns(2)
            for i, (n, img) in enumerate(items):
                with cols[i % 2]:
                    st.image(img, caption=n, width='stretch')

    if reports or datas:
        st.markdown("**Files**")
        cols = st.columns(3)
        for i, f in enumerate(reports + datas):
            with cols[i % 3]:
                with open(f, 'rb') as fh:
                    st.download_button(os.path.basename(f), fh.read(), file_name=os.path.basename(f),
                                       icon=':material/picture_as_pdf:' if f.lower().endswith('.pdf')
                                       else ':material/download:', key=f"{key}|dl|{f}", width='stretch')
    if run.get('folders'):
        st.caption("Folder: " + " · ".join(os.path.relpath(f, files.RESULTS_DIR) for f in run['folders']))


#######################################################################################################################
# PyPowerSim: helper
#######################################################################################################################
# data = {'calc', 'setup' (loaded by the toolkit), 'time', 'freq', 'sweep'} recorded by adapter.run
LOSS_PARTS = {'p_T_c': 'transistor conduction', 'p_T_s': 'transistor switching', 'p_D_c': 'diode conduction',
              'p_D_s': 'diode switching'}
TRA_LOSS = {'p_PC': 'primary copper', 'p_SC': 'secondary copper', 'p_CC': 'core'}
GROUP_TITLES = {'sw': 'switches', 'cap': 'capacitors', 'tra': 'transformer'}
# time signals of the topology classes (src/topo/class*.py out/calcTime): name -> (meaning, unit)
CHANNELS = {
    # switching (p.u. = normalised to Vdc/2 or the carrier amplitude)
    'v_a_ref': ("reference voltage phase a", 'V'), 'v_b_ref': ("reference voltage phase b", 'V'),
    'v_c_ref': ("reference voltage phase c", 'V'), 'e_a': ("back EMF phase a", 'V'), 'e_b': ("back EMF phase b", 'V'),
    'e_c': ("back EMF phase c", 'V'), 'e': ("back EMF", 'V'),
    'xA': ("modulation function a (incl. zero sequence)", 'p.u.'), 'xB': ("modulation function b", 'p.u.'),
    'xC': ("modulation function c", 'p.u.'), 'xAs': ("sampled modulation function a", 'p.u.'),
    'xBs': ("sampled modulation function b", 'p.u.'), 'xCs': ("sampled modulation function c", 'p.u.'),
    'xAsh': ("sampled and held modulation function a", 'p.u.'), 'xBsh': ("sampled and held function b", 'p.u.'),
    'xCsh': ("sampled and held function c", 'p.u.'), 'sA': ("switching function leg A (+1 upper, −1 lower)", 'p.u.'),
    'sB': ("switching function leg B", 'p.u.'), 'sC': ("switching function leg C", 'p.u.'),
    'sD': ("switching function leg D", 'p.u.'), 'n0': ("zero-sequence signal", 'p.u.'),
    'n0A': ("zero-sequence signal bridge A", 'p.u.'), 'n0B': ("zero-sequence signal bridge B", 'p.u.'),
    'c': ("carrier", 'p.u.'), 'cA': ("carrier primary bridge", 'p.u.'), 'cB': ("carrier secondary bridge", 'p.u.'),
    # inverter AC side
    'v_a0': ("leg voltage a against the DC-link midpoint", 'V'), 'v_a': ("phase voltage a (load)", 'V'),
    'v_L_a': ("voltage over the load impedance a (v_out − E)", 'V'), 'v_a_out': ("output voltage a (after filter)", 'V'),
    'v_n0': ("star-point voltage against the DC-link midpoint", 'V'), 'i_a': ("phase current a", 'A'),
    'i_b': ("phase current b", 'A'), 'i_c': ("phase current c / capacitor current", 'A'),
    'i_ref_a': ("current reference a (closed loop)", 'A'), 'i_ref_b': ("current reference b", 'A'),
    'i_ref_c': ("current reference c", 'A'),
    # inverter DC side
    'v_in': ("DC input voltage", 'V'), 'v_dc': ("DC-link capacitor voltage ripple", 'V'),
    'i_dc': ("DC-link current (bridge input)", 'A'),
    # DAB
    'v_ac_pri': ("AC voltage primary bridge", 'V'), 'v_ac_sec': ("AC voltage secondary bridge", 'V'),
    'v_ac_sec_ref': ("secondary AC voltage referred to the primary (·N)", 'V'),
    'v_L': ("voltage over the leakage inductance", 'V'), 'i_ac_pri': ("transformer current primary", 'A'),
    'i_ac_sec': ("transformer current secondary", 'A'), 'v_dc_pri': ("DC voltage primary (input)", 'V'),
    'v_dc_sec': ("DC voltage secondary (output)", 'V'), 'v_dc_pri_cap': ("input capacitor voltage ripple", 'V'),
    'v_dc_sec_cap': ("output capacitor voltage ripple", 'V'), 'i_dc_in': ("DC input current (source)", 'A'),
    'i_dc_pri': ("DC current primary bridge", 'A'), 'i_dc_sec': ("DC current secondary bridge", 'A'),
    'i_dc_out': ("DC output current (load)", 'A'), 'i_c_pri': ("input capacitor current", 'A'),
    'i_c_sec': ("output capacitor current", 'A'),
}
# spectra (src/general/calcFreq.py, calcSpec.py): name -> (meaning, unit)
SPECTRA = {'Sa': ("switching function leg A", 'p.u.'), 'Sb': ("switching function leg B", 'p.u.'),
           'Xas': ("sampled modulation function a", 'p.u.'), 'Xbs': ("sampled modulation function b", 'p.u.'),
           'I_a': ("phase current a", 'A'), 'V_a': ("phase voltage a", 'V'), 'V_a0': ("leg voltage a", 'V'),
           'I_dc': ("DC-link current", 'A'), 'V_dc': ("DC-link voltage ripple", 'V')}
# sweep quantities (src/mode/inv/calcSweep.py): suffix -> (meaning, unit factor)
SWEEP_SUFFIX = {'_v1_eff': "RMS value of the fundamental", '_eff': "RMS value",
                '_thd': "RMS value of the harmonics (distortion, without fundamental)"}


def channel(name):
    """(meaning, unit) of a time signal or spectrum."""
    if name in CHANNELS:
        return CHANNELS[name]
    if name in SPECTRA:
        return SPECTRA[name]
    return ('', 'V' if name.lower().startswith('v') else 'A' if name.lower().startswith('i') else '')


def sweep_info(q):
    """(meaning, unit) of a sweep quantity like 'V_a_v1_eff'."""
    for suf, txt in SWEEP_SUFFIX.items():
        if q.endswith(suf):
            base = q[:-len(suf)]
            what = channel(base)[0] or base
            return f"{txt} of the {what}", 'V' if base.startswith('V') else 'A'
    return '', ''


def _num(x):
    """Real scalar of a setup value (complex values as magnitude)."""
    try:
        x = complex(x)
        return abs(x) if x.imag else x.real
    except (TypeError, ValueError):
        return np.nan


def _window(t, fel, which):
    """Index range of the evaluation window: the last fundamental period (as the toolkit plots) or everything."""
    n = len(t)
    if which == 'Full simulation' or not fel or n < 3:
        return 0, n
    K = int(np.round((t[-1] - t[0]) * fel))
    if K < 2:
        return 0, n
    return int((n - 1) / K) * (K - 1), n


def _channels(time):
    """Time signals {group: {name: array}} of the switching, AC and DC side."""
    out = {}
    sw = time.get('Sw')
    if isinstance(sw, pd.DataFrame):
        out['Switching'] = {c: sw[c].to_numpy() for c in sw.columns if c != 't'}
    for g, title in (('Ac', 'AC side'), ('Dc', 'DC side')):
        if isinstance(time.get(g), dict):
            out[title] = {}
            for k, v in time[g].items():
                if isinstance(v, dict):                                # closed loop: i_ref = {'A': …, 'B': …}
                    out[title].update({f"{k}_{p.lower()}": np.asarray(x, dtype=float) for p, x in v.items()})
                else:
                    out[title][k] = np.asarray(v, dtype=float)
    return out


def _frames(group):
    """{component: DataFrame} of a loss/electrical group (the DAB transformer is a single DataFrame)."""
    if isinstance(group, pd.DataFrame):
        return {'Tr': group}
    return {k: v if isinstance(v, pd.DataFrame) else pd.DataFrame(v) for k, v in (group or {}).items()}


def _losses(time, T):
    """Mean losses per component over the window: one row per switch/capacitor/transformer."""
    rows = []
    loss = time.get('Loss')
    if not isinstance(loss, dict):
        return pd.DataFrame()
    for g, grp in loss.items():
        for comp, df in _frames(grp).items():
            df = T.cut(df)
            row = {'component': comp, 'type': GROUP_TITLES.get(g, g)}
            for c in LOSS_PARTS:
                if c in df:
                    row[LOSS_PARTS[c]] = float(df[c].mean())
            for c in df.columns:                                       # transformer: core and copper losses
                if c not in LOSS_PARTS and c not in ('p_T', 'p_D', 'p_L'):
                    row[c] = float(df[c].mean())
            row['total'] = float(df['p_L'].mean()) if 'p_L' in df else float(df.mean().sum())
            rows.append(row)
    return pd.DataFrame(rows)


def _temps(time):
    """Temperature traces {'S1 · T1': array} of all components."""
    out = {}
    ther = time.get('Ther')
    if not isinstance(ther, dict):
        return out
    for g, grp in ther.items():
        items = grp.items() if isinstance(grp, dict) else (grp.to_dict('series').items()
                                                           if isinstance(grp, pd.DataFrame) else [])
        for k, v in items:
            out[f"{GROUP_TITLES.get(g, g)} · {k}"] = np.asarray(v, dtype=float)
    return out


def _t(time, fsim):
    """
    Time axis of the AC/DC signals (most frequent length): time['Sw']['t'] if it has their length (steady state,
    transient; time['t'] has a different step there), otherwise n samples at 1/fsim (closed loop: Sw holds only the
    last control segment and time['t'] is scaled to one fundamental period although the signals cover tmax).
    """
    sw, t = time.get('Sw'), np.asarray(time.get('t'), dtype=float)
    lens = [len(v) for g in ('Ac', 'Dc') for v in (time.get(g) or {}).values() if not isinstance(v, dict)]
    n = max(set(lens), key=lens.count) if lens else len(t)
    if isinstance(sw, pd.DataFrame) and 't' in sw and len(sw) == n:
        return sw['t'].to_numpy(dtype=float)
    return np.arange(n) / fsim if fsim else t


class _Time:
    """
    Time axes and evaluation windows per array length: the arrays of one run have different lengths (transient: AC/DC
    signals of the last period only, losses and temperatures over tmax; DAB: some AC voltages finer sampled).
    """
    def __init__(self, time, fsim, fel, which):
        self.t, self.fsim, self.fel, self.which = _t(time, fsim), fsim, fel, which

    def axis(self, n):
        if n == len(self.t):
            return self.t
        return np.arange(n) / self.fsim if self.fsim else np.linspace(self.t[0], self.t[-1], n)

    def win(self, n):
        return _window(self.axis(n), self.fel, self.which)

    def cut(self, x):
        i0, i1 = self.win(len(x))
        return x.iloc[i0:i1] if isinstance(x, (pd.DataFrame, pd.Series)) else np.asarray(x)[i0:i1]


def _rms(x):
    return float(np.sqrt(np.mean(np.asarray(x, dtype=float) ** 2)))


def _power(time, s, T, p_loss):
    """
    Output power and efficiency as in the results tables of the toolkit (src/plot/plotResults*.py):
    inverter P = |Z|·I_rms²·cos(∠Z) of i_a and η = 1 − P_loss/P; DAB Pin = I_rms·V_rms at the DC input,
    Po at the DC output and η = 1 − P_loss/Pin.
    """
    ac, dc = time.get('Ac') or {}, time.get('Dc') or {}
    try:
        if s['Top']['sourceType'] == 'DAB':
            p_in = _rms(T.cut(dc['i_dc_in'])) * _rms(T.cut(dc['v_dc_pri']))
            p_out = _rms(T.cut(dc['i_dc_out'])) * _rms(T.cut(dc['v_dc_sec']))
            return p_out, 1 - p_loss / p_in
        w = 2 * np.pi * _num(s['Top']['fel'])
        R, L = _num(s['Top']['R']), _num(s['Top']['L'])
        p_out = np.hypot(R, w * L) * _rms(T.cut(ac['i_a'])) ** 2 * np.cos(np.arctan2(w * L, R))
        return p_out, 1 - p_loss / p_out
    except (KeyError, ZeroDivisionError, TypeError, ValueError):
        return np.nan, np.nan


def _thd(x, t):
    """THD of an AC signal over its window with the function of the toolkit (src/general/helpFnc.thd)."""
    try:
        from src.general.helpFnc import thd
        return float(thd(np.asarray(x, dtype=float), np.asarray(t, dtype=float), 2 / np.sqrt(2), 1))
    except Exception:
        return np.nan


#######################################################################################################################
# PyPowerSim: views
#######################################################################################################################
def _waveforms(time, T, key):
    chans = _channels(time)
    if not chans:
        st.info("No time signals in this run.")
        return
    c1, c2 = st.columns([1, 3])
    with c1:
        grp = st.selectbox("Signals", list(chans), index=min(1, len(chans) - 1), key=f"{key}|grp")
    names = list(chans[grp])
    dflt = [n for n in names if n.startswith(('i_a', 'i_b', 'i_c', 'v_a', 'i_ac', 'v_ac', 'v_dc', 'i_dc'))][:3] \
        or names[:2]
    with c2:
        pick = st.multiselect("Channels", names, default=dflt, key=f"{key}|ch|{grp}",
                              format_func=lambda n: f"{n} — {channel(n)[0]}" if channel(n)[0] else n)
    if not pick:
        return
    units = list(dict.fromkeys(channel(n)[1] for n in pick))
    right = units[1] if len(units) > 1 else None                      # second unit on the right axis
    fig = go.Figure()
    rows = []
    for i, n in enumerate(pick):
        x = chans[grp][n]
        tt = T.axis(len(x))
        a, b = T.win(len(x))
        s = downsample(b - a)
        u = channel(n)[1]
        fig.add_trace(line(tt[a:b:s] * 1e3, x[a:b:s], f"{n} [{u}]", i, width=1.5,
                           yaxis='y2' if right and u == right else 'y'))
        w = x[a:b]
        rows.append({'signal': n, 'meaning': channel(n)[0], 'unit': u, 'mean': float(np.mean(w)),
                     'RMS': float(np.sqrt(np.mean(w ** 2))), 'min': float(np.min(w)), 'max': float(np.max(w))})
    fig.update_layout(xaxis_title="time t [ms]", yaxis_title=f"[{units[0]}]", hovermode='x unified')
    if right:
        fig.update_layout(yaxis2=dict(title=f"[{right}]", overlaying='y', side='right', showgrid=False))
    if len(units) > 2:
        st.caption(f":material/warning: More than two units selected ({', '.join(units)}); all but '{right}' share "
                   f"the left axis.")
    show(fig, height=420, key=f"{key}|wave")
    table(pd.DataFrame(rows), hide_index=True,
          column_config={c: st.column_config.NumberColumn(c, format='%.4g', help=f"{c} over the evaluation window, "
                                                                                  f"in the unit of the signal")
                         for c in ('mean', 'RMS', 'min', 'max')})


def _loss_view(time, T, key):
    df = _losses(time, T)
    if df.empty:
        st.info("No losses in this run (Par › PWM › Losses = 0 gives ideal switches).")
        return
    parts = [c for c in df.columns if c not in ('component', 'type', 'total')]
    fig = go.Figure([go.Bar(x=df['component'], y=df[c], name=TRA_LOSS.get(c, c), marker_color=color(i))
                     for i, c in enumerate(parts)])
    fig.update_layout(barmode='stack', yaxis_title="mean loss P [W]", xaxis_title=None)
    show(fig, height=380, key=f"{key}|bar")
    out = df.rename(columns={**{c: f"{TRA_LOSS.get(c, c)} [W]" for c in parts}, 'total': 'total [W]'})
    out.insert(len(out.columns), 'share [%]', 100 * df['total'] / df['total'].sum())
    table(out, hide_index=True)
    st.caption("Mean values over the evaluation window. Switches S1 … Sn: conduction losses (on-state voltage × "
               "current, from the Vf/Vfd tables or Ron) and switching losses (Eon, Eoff, Erec scaled to the actual "
               "voltage and current, × switching events) of transistor and anti-parallel diode. Capacitors: ESR losses "
               "of the ripple current. Transformer: copper losses of both windings and core losses.")


def _temp_view(time, T, key):
    tr = {n: T.cut(v) for n, v in _temps(time).items()}
    if not tr:
        st.info("No temperatures in this run (no losses, or not stored by this mode).")
        return
    names = list(tr)
    dflt = [n for n in names if n.endswith(('· T1', '· D1'))] or names[:2]
    pick = st.multiselect("Components", names, default=dflt, key=f"{key}|th")
    if pick:
        fig = go.Figure()
        for i, n in enumerate(pick):
            a, b = T.win(len(_temps(time)[n]))
            tt = T.axis(len(_temps(time)[n]))[a:b]
            s = downsample(len(tt))
            fig.add_trace(line(tt[::s] * 1e3, tr[n][::s], n, i, width=1.5))
        fig.update_layout(xaxis_title="time [ms]", yaxis_title="temperature [°C]", hovermode='x unified')
        show(fig, height=380, key=f"{key}|temp")
    table(pd.DataFrame([{'component': n, 'meaning': _temp_meaning(n), 'mean [°C]': float(np.mean(v)),
                         'max [°C]': float(np.max(v)), 'swing ΔT [K]': float(np.max(v) - np.min(v))}
                        for n, v in tr.items()]), hide_index=True)
    st.caption("Values over the evaluation window. The swing ΔT = max − min is the temperature cycle that drives the "
               "thermo-mechanical fatigue of the devices (solder, bond wires).")


def _temp_meaning(n):
    k = n.split('· ')[-1]
    if k[:1] == 'T' and k[1:].isdigit():
        return f"junction temperature transistor {k[1:]}"
    if k[:1] == 'D' and k[1:].isdigit():
        return f"junction temperature diode {k[1:]}"
    if k[:1] == 'C' and k[1:].isdigit():
        return f"case temperature {'capacitor' if n.startswith('capacitors') else 'switch'} {k[1:]}"
    return {'PC': 'primary winding', 'SC': 'secondary winding', 'CC': 'core'}.get(k, '')


def _spectrum(freq, setup, key):
    groups = {g: v for g, v in (freq or {}).items() if isinstance(v, dict) and v}
    if not groups:
        st.info("No spectra in this run.")
        return
    fsim, fel = _num(setup['Exp']['fsim']), _num(setup['Top']['fel'])
    c1, c2, c3 = st.columns([1, 2, 1])
    with c1:
        g = st.selectbox("Side", list(groups), key=f"{key}|fg")
    with c2:
        n = st.selectbox("Signal", list(groups[g]), key=f"{key}|fs|{g}",
                         format_func=lambda x: f"{x} — {channel(x)[0]}" if channel(x)[0] else x)
    with c3:
        H = st.number_input("Harmonics", 10, 2000, 50, step=10, key=f"{key}|fh",
                            help="Highest harmonic order shown")
    u = channel(n)[1]
    y = np.abs(np.asarray(groups[g][n], dtype=float))
    h = fsim * np.linspace(0, 0.5, len(y)) / fel                    # harmonic order, as in the toolkit plots
    m = h <= H
    log = st.toggle("Logarithmic amplitude", value=True, key=f"{key}|flog")
    ys = y[m]
    base = max(np.nanmax(ys) * 1e-6, np.nanmin(ys[ys > 0])) if log and np.any(ys > 0) else 0.0
    keep = ys >= base if log else ys > 0                              # skip empty bins (dots on the axis)
    fig = go.Figure(stem(h[m][keep], ys[keep], n, 0, f"order %{{x:.2f}}<br>%{{y:.4g}} {u}", base=base))
    fig.update_layout(xaxis_title="harmonic order h = f / fel", yaxis_title=f"amplitude [{u}]", showlegend=False)
    if log:
        fig.update_yaxes(type='log')
    show(fig, height=380, key=f"{key}|spec")
    # largest components
    idx = np.where(m & (h > 0.5))[0]
    top = idx[np.argsort(y[idx])[::-1][:10]] if idx.size else []
    one = y[np.argmin(np.abs(h - 1))] if np.any(np.abs(h - 1) < 0.5) else np.nan
    rows = [{'order h': float(h[i]), 'frequency [Hz]': float(h[i] * fel), f'amplitude [{u}]': float(y[i]),
             'of fundamental [%]': 100 * float(y[i] / one) if one else np.nan} for i in sorted(top)]
    if rows:
        st.caption(f"Ten largest components up to order {H} (fundamental fel = {fmt(fel, 'Hz')}). Amplitudes as "
                   f"computed by the toolkit (src/general/calcFreq.py, peak values of the single-sided spectrum).")
        table(pd.DataFrame(rows), hide_index=True)


DERIVED = {'Po': ("output active power", 'W'), 'Qo': ("output reactive power", 'VAr'),
           'So': ("output apparent power", 'VA'), 'Io': ("output current (RMS, complex)", 'A'),
           'Vo': ("output voltage (RMS, complex)", 'V'), 'PhiVI': ("angle between voltage and current", 'rad'),
           'Mi': ("modulation index", '–'), 'PhiDAB': ("phase shift of the DAB", '°')}


def _operating_point(setup):
    from gui.core import forms
    names = {}
    for sec in ('Dat', 'Top', 'Exp'):
        for p in forms.schema().get(sec, []):
            names[(sec, p['var'])] = (p['name'] or p['var'], p['unit'] or '')
    rows = []
    for sec, cat in (('Dat', 'stat'), ('Top', None), ('Exp', None)):
        d = setup[sec][cat] if cat else setup[sec]
        for k, v in d.items():
            if isinstance(v, (int, float, complex, np.number)) and not isinstance(v, bool):
                c = complex(v)
                meaning, unit = DERIVED.get(k) if k in DERIVED else names.get((sec, k), ('', ''))
                rows.append({'parameter': f"{sec}.{cat + '.' if cat else ''}{k}", 'meaning': meaning,
                             'value': abs(c) if c.imag else c.real, 'unit': unit,
                             'angle [°]': float(np.degrees(np.angle(c))) if c.imag else None})
    st.caption("Values of the setup as used by the simulation, after the toolkit derived the operating point "
               "(src/general/genLoadInput.py: Po, Qo, So, Io, Vo, PhiVI from the controlled quantity). Complex "
               "quantities are shown as magnitude and angle.")
    table(pd.DataFrame(rows), hide_index=True, height=420,
          column_config={'value': st.column_config.NumberColumn(format='%.5g')})


def pps_run(data, setup, key):
    """Steady-state, transient and closed-loop results: key figures, waveforms, losses, temperatures, spectra."""
    time, s = data.get('time') or {}, data.get('setup') or {}
    fel = _num(s['Top']['fel'])
    steady = data.get('calc', '').startswith(('calcSteady', 'calcSweep'))
    which = st.segmented_control("Evaluation window", ['Last period', 'Full simulation'],
                                 default='Last period' if steady else 'Full simulation', key=f"{key}|win") or \
        'Last period'
    T = _Time(time, _num(s['Exp']['fsim']), fel, which)
    t = T.t
    i0, i1 = T.win(len(t))
    if data.get('calc', '').startswith('calcTrans'):
        st.caption("Transient: the waveforms and spectra cover the last fundamental period, losses and temperatures "
                   "the whole simulation time.")

    # Key figures
    loss = _losses(time, T)
    p_loss = float(loss['total'].sum()) if not loss.empty else np.nan
    p_out, eff = _power(time, s, T, p_loss)
    t_max = max((float(np.max(T.cut(v))) for v in _temps(time).values()), default=np.nan)
    ac = time.get('Ac') or {}
    dab = s['Top']['sourceType'] == 'DAB'
    c = st.columns(5)
    c[0].metric("Output power", fmt(p_out, 'W'),
                help="DAB: I_rms·V_rms at the DC output. Inverter: |Z|·I_rms²·cos∠Z of phase a (as in the results "
                     "table of the toolkit).")
    c[1].metric("Total losses", fmt(p_loss, 'W'), help="Sum of the mean losses of all components over the window.")
    c[2].metric("Efficiency", fmt(100 * eff, '%', 4) if np.isfinite(eff) else '–',
                help="1 − P_loss / P (DAB: P = DC input power), as in the results table of the toolkit.")
    c[3].metric("Max. temperature", fmt(t_max, '°C'), help="Highest temperature of all components in the window.")
    one_period = abs((t[i1 - 1] - t[i0]) * fel - 1) < 0.05 if i1 - i0 > 2 and fel else False
    if not dab and one_period and 'i_a' in ac and len(ac['i_a']) == len(t):
        c[4].metric("THD i_a", fmt(100 * _thd(np.asarray(ac['i_a'])[i0:i1], t[i0:i1]), '%'),
                    help="Total harmonic distortion of the phase current over the window (src/general/helpFnc.thd, "
                         "fundamental = first harmonic of the window, so use 'Last period').")

    tabs = st.tabs(["Waveforms", "Losses", "Temperatures", "Spectrum", "Operating point"])
    with tabs[0]:
        _waveforms(time, T, key)
    with tabs[1]:
        _loss_view(time, T, key)
    with tabs[2]:
        _temp_view(time, T, key)
    with tabs[3]:
        _spectrum(data.get('freq'), s, key)
    with tabs[4]:
        _operating_point(s)


def pps_sweep(data, setup, key):
    """Sweep results: RMS values, fundamentals and THD over the swept parameter, numerical vs analytical."""
    sw, s = data.get('sweep') or {}, data.get('setup') or {}
    xname = next((k for k in ('Mi', 'Phi', 'PhiDAB', 'phi') if k in sw), None)
    groups = {g: v for g, v in sw.items() if isinstance(v, dict) and 'num' in v}
    if xname is None or not groups:
        st.info("No sweep data in this run.")
        generic(data, setup, key)
        return
    x = np.asarray(sw[xname], dtype=float)
    tabs = st.tabs(["Sweep", "Waveforms (last point)", "Spectrum (last point)"])
    with tabs[0]:
        c1, c2 = st.columns([1, 3])
        with c1:
            g = st.selectbox("Side", list(groups), key=f"{key}|sg")
        cols = list(groups[g]['num'].columns)
        with c2:
            pick = st.multiselect("Quantities", cols, default=[c for c in cols if c.endswith('_thd')][:2] or cols[:2],
                                  key=f"{key}|sq|{g}", format_func=lambda q: f"{q} — {sweep_info(q)[0]}")
        unit = lambda q: sweep_info(q)[1]
        scale = lambda q: 1.0                                          # '_thd' of the toolkit is an RMS value [V, A]
        units = list(dict.fromkeys(unit(q) for q in pick))
        right = units[1] if len(units) > 1 else None
        fig = go.Figure()
        for i, q in enumerate(pick):
            ax = 'y2' if right and unit(q) == right else 'y'
            fig.add_trace(line(x, scale(q) * groups[g]['num'][q], f"{q} [{unit(q)}] simulation", i,
                               mode='lines+markers', yaxis=ax))
            if 'ana' in groups[g] and q in groups[g]['ana']:
                fig.add_trace(line(x, scale(q) * groups[g]['ana'][q], f"{q} [{unit(q)}] analytical", i, dash='dash',
                                   yaxis=ax))
        xlab = {'Mi': "modulation index Mi [–]", 'PhiDAB': "phase shift [°]"}.get(xname, xname)
        fig.update_layout(xaxis_title=xlab, yaxis_title=f"[{units[0]}]" if units else None, hovermode='x unified')
        if right:
            fig.update_layout(yaxis2=dict(title=f"[{right}]", overlaying='y', side='right', showgrid=False))
        show(fig, height=420, key=f"{key}|sweep")
        st.caption("Solid: numerical simulation, dashed: analytical solution of the toolkit.")
        tab = pd.DataFrame({xlab: x})
        for q in cols:
            tab[f"{q} [{unit(q)}]"] = scale(q) * groups[g]['num'][q].to_numpy(dtype=float)
        table(tab, hide_index=True)
        st.caption("Columns: " + "; ".join(f"**{q}** {sweep_info(q)[0]} [{unit(q)}]" for q in cols) +
                   ". Simulated values over the swept operating points (src/general/calcDistNum.py); the distortion "
                   "'_thd' is the RMS value of all harmonics √(X_rms² − X_1,rms²), the ratio THD = X_thd / X_1,rms.")
    time = data.get('time') or {}
    with tabs[1]:
        _waveforms(time, _Time(time, _num(s['Exp']['fsim']), _num(s['Top']['fel']), 'Last period'), f"{key}|w")
    with tabs[2]:
        _spectrum(data.get('freq'), s, f"{key}|f")


#######################################################################################################################
# Interactive results per mode: mode -> function(data, setup, key)
#######################################################################################################################
RENDER = {0: pps_sweep, 1: pps_run, 2: pps_run, 3: pps_run}


#######################################################################################################################
# Run
#######################################################################################################################
def run_view(run, key='res'):
    """Interactive results (if the data are available), figures and files, and the log of one run."""
    if run is None:
        st.info("No results yet — start a simulation in the **Run** view or open a results folder.")
        return
    if run.get('error'):
        st.error(f"The simulation stopped with an error:\n\n```\n{run['error']}\n```")
    st.caption(f"Run {run.get('time', '')} · {fmt(run.get('duration', 0))} s")
    view = st.segmented_control("Show", ['Interactive', 'Figures & files', 'Console log'], default='Interactive',
                                key=f"{key}|v", label_visibility='collapsed')
    if view == 'Interactive':
        if run.get('data') is not None:
            try:
                RENDER.get(run['mode'], generic)(run['data'], run.get('setup'), key=f"{key}|{run['mode']}")
            except Exception as e:
                st.warning(f"Interactive charts not available for this run ({type(e).__name__}: {e}). The figures "
                           f"of the toolkit are shown under 'Figures & files'.")
        else:
            st.info("No interactive data stored for this run — see 'Figures & files'.")
    elif view == 'Figures & files':
        gallery(run, key=f"{key}|g")
    elif view == 'Console log':
        st.code(run.get('log') or '(log not stored for this run)', language=None, height=500)
