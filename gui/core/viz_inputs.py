#######################################################################################################################
#######################################################################################################################
# Title:        Power Electronics Simulation Toolkit (PyPowerSim)
# Topic:        Power Electronics
# File:         viz_inputs.py
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
This function shows the input workbooks. The device parameter workbooks (para/Swi, para/Cap, para/Tra) are shown with
dedicated views — parameter table, characteristics (the look-up tables of the loss and component models as charts),
thermal Foster networks (Zth curves) — and can be edited and saved as in the RAISE GUI. Other workbooks use the
generic viewer (one tab per sheet: table, and a chart of the numeric columns over the first column).
INPUTS (mode -> function(cfg)) holds dedicated input views per operating mode for the Inputs view of a mode page.
"""

#######################################################################################################################
# Import libs
#######################################################################################################################
from gui import adapter as A
from gui.core import devices as D
from gui.core import files, forms
from gui.core.style import show, line, downsample, table, fmt, color, num_format

import copy
import os

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st


#######################################################################################################################
# Generic viewer
#######################################################################################################################
def arrow_safe(df):
    """Mixed text/number columns (typical for parameter sheets) as text, so the table renders without warnings."""
    df = df.copy()
    df.columns = [str(c) for c in df.columns]
    for c in df.columns:
        if df[c].dtype == object and not all(isinstance(v, str) for v in df[c].dropna()):
            df[c] = df[c].map(lambda v: '' if v is None or (isinstance(v, float) and v != v) else str(v))
    return df


def sheet(df, key):
    """Table of a sheet and a line chart of its numeric columns over the first column."""
    df = arrow_safe(df.dropna(how='all').dropna(axis=1, how='all'))
    if df.empty:
        st.info("Empty sheet.")
        return
    num = [c for c in df.columns[1:] if pd.api.types.is_numeric_dtype(df[c])]
    x = df.columns[0]
    if num and len(df) > 5 and pd.api.types.is_numeric_dtype(df[x]):
        pick = st.multiselect("Columns", num, default=num[:min(4, len(num))], key=f"{key}|cols")
        if pick:
            s = downsample(len(df))
            fig = go.Figure([line(df[x].iloc[::s], df[c].iloc[::s], str(c), i) for i, c in enumerate(pick)])
            fig.update_layout(xaxis_title=str(x))
            show(fig, key=f"{key}|ts")
    with st.expander(f"Table ({len(df)} rows × {len(df.columns)} columns)", expanded=not num):
        table(df, hide_index=True, height=320)


def workbook(kind, name, key, field=None):
    """An input workbook: device view for para/Swi, para/Cap, para/Tra, otherwise one tab per sheet."""
    if not name or name not in files.list_files(kind):
        st.info(f"Workbook '{name}' not found in {A.FOLDERS.get(kind, kind)}/.")
        return
    if D.kind_of(kind):
        device(kind, name, key, field)
        return
    try:
        sheets = files.load_workbook(kind, name, files.mtime(kind, name))
    except Exception as e:
        st.error(f"The workbook could not be read: {e}")
        return
    names = [s for s in sheets if s.lower() != 'readme'] or list(sheets)
    for tab, s in zip(st.tabs(names), names):
        with tab:
            sheet(sheets[s], key=f"{key}|{s}")


#######################################################################################################################
# Device workbooks: helper
#######################################################################################################################
EDIT_COLS = ['Typical'] + D.VALUES
SI = '~s'                                                             # Plotly SI prefixes (5m = 5e-3)
MODEL_TEXT = {'Constant': 'constant value (Typical)', 'Linear': 'linear model (Typical)',
              'Tabular': 'look-up table (vector Value-1 … 10 / LUT)', 'General': 'general rating',
              'Losses': 'loss model', '1D-Foster': 'Foster network branches (Value-1 … 10)'}


def _base(kind, name):
    """Workbook as read from disk, cached in the session per modification time."""
    k = f"dev|{kind}|{name}|{files.mtime(kind, name)}"
    if k not in st.session_state:
        st.session_state[k] = D.read(kind, name)
    return st.session_state[k]


def _ver(kind, name):
    return st.session_state.get(f"devver|{kind}|{name}", 0)


def _reset(kind, name):
    st.session_state[f"devver|{kind}|{name}"] = _ver(kind, name) + 1
    for k in [k for k in st.session_state if str(k).startswith(f"dev|{kind}|{name}|")]:
        del st.session_state[k]


def _vec(elec, sym):
    row = elec[elec['Symbol'] == sym]
    if row.empty:
        return []
    return [v for v in pd.to_numeric(row.iloc[0][D.VALUES], errors='coerce').tolist() if np.isfinite(v)]


def _values(row):
    v = [x for x in pd.to_numeric(row[D.VALUES], errors='coerce').tolist() if np.isfinite(x)]
    return ', '.join(fmt(x, digits=4) for x in v)


def _numeric(df, cols):
    df = df.copy()
    for c in cols:
        df[c] = pd.to_numeric(df[c], errors='coerce')
    return df


#######################################################################################################################
# Device workbooks: views
#######################################################################################################################
def _header(dev):
    rd = dev['readme']
    tags = [str(v) for k, v, _ in rd if k in ('Type', 'PNR', 'Package', 'Qualification') and v not in (None, '')]
    st.markdown(' '.join(f":blue-badge[{t}]" for t in tags) or ' ')
    ratings = [(k, v, u) for k, v, u in rd if k not in ('Type', 'PNR', 'Package', 'Qualification')]
    if ratings:
        cols = st.columns(len(ratings) + 1)
        for c, (k, v, u) in zip(cols, ratings):
            c.metric(k, fmt(v, u or '') if isinstance(v, (int, float)) else (str(v) if v not in (None, '') else '–'))
    info = dev['info']
    src = info.get('Datasheet')
    st.caption(f"Author: {info.get('Author', '–')} · date: {str(info.get('Date', '–'))[:10]}" +
               (f" · [datasheet]({src})" if isinstance(src, str) and src.startswith('http') else
                (f" · datasheet: {src}" if src else '')))


def _param_table(df, title):
    """Parameter table with value, unit, model, vector and description."""
    rows = []
    for _, r in df.iterrows():
        if r['Symbol'] in (None, '') or (isinstance(r['Symbol'], float) and np.isnan(r['Symbol'])):
            continue
        typ = pd.to_numeric(pd.Series([r['Typical']]), errors='coerce').iloc[0]
        rows.append({'Symbol': r['Symbol'], 'Parameter': r['Parameter'], 'Typical': typ,
                     'Unit': r['Unit'] if isinstance(r['Unit'], str) else '',
                     'Model': MODEL_TEXT.get(r['Model'], r['Model']) if isinstance(r['Model'], str) else '',
                     'Vector (Value-1 … 10)': _values(r),
                     'Description': r['Description'] if isinstance(r['Description'], str) else ''})
    st.markdown(f"**{title}**")
    table(pd.DataFrame(rows), hide_index=True, width='stretch',
          column_config={'Typical': st.column_config.NumberColumn("Typical value", format='%.4g',
                                                                  help="Value of the constant/linear models; "
                                                                       "for vectors the typical operating point"),
                         'Vector (Value-1 … 10)': st.column_config.TextColumn(
                             help="Axis points of the look-up tables (e.g. If, Tj) or Foster branches"),
                         'Description': st.column_config.TextColumn(width='large')})


def _lut_chart(sym, lut, key):
    df = D.lut_frame(lut)
    fig = go.Figure()
    for i, c in enumerate(df.columns):
        fig.add_trace(line(df.index, df[c], c, i, mode='lines+markers'))
    fig.update_layout(title=dict(text=lut['title'], font=dict(size=13)), xaxis_title=lut['xlabel'],
                      yaxis_title=f"{sym} [{lut['unit']}]", hovermode='x unified')
    fig.update_yaxes(tickformat=SI)
    show(fig, height=320, key=key, legend_below=True)


def device_luts(dev, key):
    """Characteristics: every defined look-up table as curves over its x-axis, one curve per temperature."""
    luts = {s: l for s, l in dev['luts'].items() if np.isfinite(D.lut_frame(l).to_numpy(dtype=float)).any()}
    empty = [s for s in dev['luts'] if s not in luts]
    if not luts:
        st.info("No look-up tables are defined in this workbook.")
        return
    st.caption("Look-up tables of the electrical sheet as read by the loader of the toolkit: rows = values of the "
               "x-axis vector (parameter table, e.g. If), columns = junction temperatures (vector Tj). Between the "
               "points the toolkit interpolates linearly (scipy griddata / RegularGridInterpolator).")
    syms = list(luts)
    for i in range(0, len(syms), 2):
        cols = st.columns(2)
        for c, s in zip(cols, syms[i:i + 2]):
            with c:
                _lut_chart(s, luts[s], key=f"{key}|lut|{s}")
    if empty:
        st.caption("Not defined (empty tables): " + ', '.join(empty) + "." +
                   (" The capacitances are needed only for MOSFETs with switching losses from capacitances "
                    "(Par › PWM › swloss = 1)." if dev['kind'] == 'swi' else ''))
    with st.expander("Values of the look-up tables", icon=':material/table:'):
        s = st.selectbox("Table", syms, key=f"{key}|lutsel",
                         format_func=lambda x: f"{x} — {luts[x]['title']}")
        df = D.lut_frame(luts[s])
        df.index.name = luts[s]['xlabel']
        st.caption(f"Values in {luts[s]['unit']}; rows: {luts[s]['xlabel']}, columns: junction temperature.")
        table(df.reset_index(), hide_index=True)


def device_thermal(dev, key):
    """Thermal Foster networks: Zth(t) of every path and the branches with time constants."""
    net = {p: rc for p, rc in D.foster(dev).items() if rc[0]}
    if not net:
        st.info("No Foster network defined.")
        return
    st.caption("Foster networks of the thermal sheet (branches i = Value-1 … 10). The toolkit adds the case → "
               "ambient (heatsink) network according to Par › Ther › Heatsink. Step response for a constant loss "
               "P: ΔT(t) = P · Zth(t).")
    st.latex(r"Z_{th}(t) = \sum_i R_{th,i}\left(1 - e^{-t/\tau_i}\right), \qquad \tau_i = R_{th,i}\,C_{th,i}")
    tau = [r * c for R, C in net.values() for r, c in zip(R, C) if r > 0 and c > 0]
    t = np.logspace(np.floor(np.log10(min(tau) / 10)) if tau else -6, np.ceil(np.log10(max(tau) * 10)) if tau else 3,
                    300)
    fig = go.Figure()
    for i, (p, (R, C)) in enumerate(net.items()):
        fig.add_trace(line(t, D.zth(R, C, t), f"{p}: {D.THER_PATHS.get(p, p)} (Rth = {fmt(sum(R), 'K/W')})", i))
    fig.update_layout(xaxis_title="time t [s]", yaxis_title="thermal impedance Zth [K/W]", hovermode='x unified')
    fig.update_xaxes(type='log', tickformat=SI)
    show(fig, height=360, key=f"{key}|zth", legend_below=True)
    rows = [{'path': f"{p} ({D.THER_PATHS.get(p, p)})", 'branch': i + 1, 'Rth [K/W]': r, 'Cth [Ws/K]': c,
             'τ = Rth·Cth [s]': r * c} for p, (R, C) in net.items() for i, (r, c) in enumerate(zip(R, C))]
    table(pd.DataFrame(rows), hide_index=True,
          column_config={'Rth [K/W]': st.column_config.NumberColumn(format='%.4g'),
                         'Cth [Ws/K]': st.column_config.NumberColumn(format='%.4g'),
                         'τ = Rth·Cth [s]': st.column_config.NumberColumn(format='%.3g')})
    neg = [p for p, (R, C) in net.items() if min(R + C) < 0]
    if neg:
        st.caption(f":material/warning: Negative branch values in {', '.join(neg)} — a fitted network can have "
                   f"negative elements; only the total Zth(t) is physical.")


#######################################################################################################################
# Device workbooks: editor
#######################################################################################################################
def _editor(base, kind, name, key):
    """Editors of the parameter tables, LUTs and Foster networks; returns the edited device dict."""
    v = _ver(kind, name)
    cur = copy.deepcopy(base)
    st.caption("Edit the values in the tables (double-click a cell). The charts of the other tabs show the edited "
               "values immediately; the workbook is changed only with **Save workbook**. Axis values of the look-up "
               "tables (e.g. If, Tj) are the vectors Value-1 … 10 of the parameter table.")
    t1, t2, t3 = st.tabs(["Parameter table", "Look-up tables", "Thermal network"])
    cc = {c: st.column_config.NumberColumn(c, format='%.4g') for c in EDIT_COLS}
    cc.update({'Symbol': st.column_config.TextColumn(disabled=True), 'Unit': st.column_config.TextColumn(disabled=True),
               'Parameter': st.column_config.TextColumn(disabled=True, width='medium')})
    with t1:
        df = _numeric(base['elec'], EDIT_COLS)[['Symbol', 'Parameter', 'Unit'] + EDIT_COLS]
        ed = st.data_editor(df, hide_index=True, column_config=cc, key=f"{key}|ed|elec|{v}", height=380)
        cur['elec'][EDIT_COLS] = ed[EDIT_COLS].to_numpy()
    with t2:
        syms = list(base['luts'])
        x_of = {s: _vec(cur['elec'], l['xsym'])[:D.N_LUT] for s, l in base['luts'].items()}
        tj = _vec(cur['elec'], 'Tj')[:D.N_LUT]
        for tab, s in zip(st.tabs(syms), syms):
            with tab:
                lut = base['luts'][s]
                x = x_of[s]
                if not x or not tj:
                    st.info(f"Define the vectors {lut['xsym']} and Tj in the parameter table first.")
                    continue
                z = np.array([[np.nan if lut['z'][i][j] is None else lut['z'][i][j] for j in range(len(tj))]
                              for i in range(len(x))], dtype=float)
                df = pd.DataFrame(z, columns=[f"Tj = {t:g} °C" for t in tj])
                df.insert(0, lut['xlabel'], x)
                st.caption(f"{lut['title']} — values in {lut['unit']}")
                ed = st.data_editor(df, hide_index=True, key=f"{key}|ed|lut|{s}|{len(x)}x{len(tj)}|{v}",
                                    column_config={lut['xlabel']: st.column_config.NumberColumn(disabled=True,
                                                                                                format='%g'),
                                                   **{c: st.column_config.NumberColumn(format='%.4g')
                                                      for c in df.columns[1:]}})
                zz = [list(r) for r in lut['z']]
                for i in range(len(x)):
                    for j in range(len(tj)):
                        val = ed.iloc[i, j + 1]
                        zz[i][j] = None if pd.isna(val) else float(val)
                cur['luts'][s]['z'] = zz
    with t3:
        df = _numeric(base['ther'], EDIT_COLS)[['Symbol', 'Parameter', 'Unit'] + EDIT_COLS]
        ed = st.data_editor(df, hide_index=True, column_config=cc, key=f"{key}|ed|ther|{v}")
        cur['ther'][EDIT_COLS] = ed[EDIT_COLS].to_numpy()
        st.caption("Rth in K/W, Cth in Ws/K; Typical holds the total of the network.")
    return cur


def _changed(base, cur):
    a = _numeric(base['elec'], EDIT_COLS)[EDIT_COLS].to_numpy(dtype=float)
    b = _numeric(cur['elec'], EDIT_COLS)[EDIT_COLS].to_numpy(dtype=float)
    c = _numeric(base['ther'], EDIT_COLS)[EDIT_COLS].to_numpy(dtype=float)
    d = _numeric(cur['ther'], EDIT_COLS)[EDIT_COLS].to_numpy(dtype=float)
    z = any(np.array(base['luts'][s]['z'], dtype=float).tobytes() != np.array(cur['luts'][s]['z'], dtype=float)
            .tobytes() for s in base['luts'])
    return not (np.array_equal(a, b, equal_nan=True) and np.array_equal(c, d, equal_nan=True)) or z


def _save_box(kind, name, cur, changed, key, field):
    cfg = forms.cfg()
    c1, c2, c3, c4 = st.columns([3, 2, 2, 2], vertical_alignment='bottom')
    with c1:
        dst = st.text_input("Save as", value=f"{name}_edit", key=f"{key}|dst",
                            help=f"Workbook name in {A.FOLDERS[kind]}/ (without .xlsx). Use the original name to "
                                 f"overwrite it.")
    dst = ''.join(ch if ch.isalnum() or ch in '-_.' else '_' for ch in dst.strip())
    exists = dst in files.list_files(kind)
    with c2:
        ok = st.checkbox("Overwrite existing", key=f"{key}|ow", disabled=not exists,
                         help="The workbook exists already — confirm to replace it.")
    with c3:
        save = st.button("Save workbook", type='primary', icon=':material/save:', key=f"{key}|save",
                         disabled=not dst or (exists and not ok), width='stretch')
    with c4:
        if st.button("Discard changes", icon=':material/undo:', key=f"{key}|discard", disabled=not changed,
                     width='stretch'):
            _reset(kind, name)
            st.rerun()
    st.caption(":material/edit: unsaved changes" if changed else "No changes yet.")
    st.caption("Saved workbooks contain values: formulas of the source workbook (e.g. LUT axes '=F$16' or scaled "
               "values) are replaced by their current values, because the loaders of the toolkit read cached cell "
               "values. The saved workbook is checked with the loader of the toolkit.")
    if save:
        out, para, log = D.save(cur, dst, cfg)
        st.session_state[f"devsaved|{kind}"] = (dst, para is not None, log)
        if dst == name:
            _reset(kind, name)
        st.rerun()
    res = st.session_state.get(f"devsaved|{kind}")
    if res:
        dst, ok, log = res
        if ok:
            st.success(f"Saved {A.FOLDERS[kind]}/{dst}.xlsx and loaded it with the loader of the toolkit.",
                       icon=':material/check_circle:')
        else:
            st.error(f"Saved {A.FOLDERS[kind]}/{dst}.xlsx, but the loader of the toolkit failed — see the log.")
        with st.expander("Loader output"):
            st.code(log or '(no output)', language=None)
        fields = [f for f in A.RUN_FIELDS if f[1] == f"file:{kind}"]
        if fields and dst in files.list_files(kind):
            c1, c2 = st.columns([3, 2], vertical_alignment='bottom')
            with c1:
                use = st.multiselect(f"Use '{dst}' for", [f[0] for f in fields], default=[field] if field in
                                     [f[0] for f in fields] else [], key=f"{key}|use",
                                     format_func=lambda x: next(f[2] for f in fields if f[0] == x))
            with c2:
                if st.button("Use in the configuration", icon=':material/input:', key=f"{key}|usebtn",
                             disabled=not use, width='stretch'):
                    for var in use:
                        cfg[A.RUN_SECTION][var] = dst
                        st.session_state.pop(f"g|{var}", None)
                    st.session_state.pop(f"devsaved|{kind}", None)
                    st.toast(f"Configuration uses {dst}", icon=':material/check:')
                    st.rerun()


def device(kind, name, key, field=None):
    """Device workbook: header, parameters, characteristics, thermal network and the editor."""
    try:
        base = _base(kind, name)
    except Exception as e:
        st.error(f"The workbook could not be read: {e}")
        return
    st.markdown(f"#### {D.LAYOUT[kind]['title']} · {name}")
    _header(base)
    tabs = st.tabs([":material/list: Parameters", ":material/show_chart: Characteristics",
                    ":material/thermostat: Thermal network", ":material/edit: Edit & save"])
    with tabs[3]:                                                     # first: the other tabs show the edited values
        cur = _editor(base, kind, name, key)
        changed = _changed(base, cur)
        _save_box(kind, name, cur, changed, key, field)
    note = ":orange-badge[:material/edit: showing unsaved edits]" if changed else None
    with tabs[0]:
        if note:
            st.markdown(note)
        _param_table(cur['elec'], "Electrical parameters")
        _param_table(cur['ther'], "Thermal parameters")
    with tabs[1]:
        if note:
            st.markdown(note)
        device_luts(cur, key)
    with tabs[2]:
        if note:
            st.markdown(note)
        device_thermal(cur, key)


#######################################################################################################################
# Inputs per mode: mode -> function(cfg)
#######################################################################################################################
INPUTS = {}
