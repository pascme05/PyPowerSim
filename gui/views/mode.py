#######################################################################################################################
#######################################################################################################################
# Title:        Power Electronics Simulation Toolkit (PyPowerSim)
# Topic:        Power Electronics
# File:         mode.py
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
Page of one operating mode with the views workflow (if defined), parameters, inputs, run and results.
"""

#######################################################################################################################
# Import libs
#######################################################################################################################
from gui import adapter as A
from gui.core import files, forms, theory, workflow
from gui.core import viz_inputs as vi
from gui.core import viz_results as vr
from gui.views.project import files_form, file_groups, active_files

import os

import streamlit as st

VIEW_ICONS = {'Workflow': ':material/account_tree:', 'Parameters': ':material/tune:', 'Inputs': ':material/table_chart:',
              'Theory': ':material/function:', 'Run': ':material/play_arrow:', 'Results': ':material/monitoring:'}


def views(m):
    return (['Workflow'] if m in workflow.WORKFLOWS else []) + ['Parameters', 'Inputs'] + \
        (['Theory'] if m in theory.THEORY else []) + ['Run', 'Results']


#######################################################################################################################
# Theory
#######################################################################################################################
def theory_view(m):
    theory.THEORY[m](forms.cfg())


#######################################################################################################################
# Parameters
#######################################################################################################################
def parameters(m):
    d = A.MODES[m]
    c1, c2 = st.columns([5, 1])
    with c2:
        forms.reset_section(d['section'])
    forms.section(d['section'])
    for sec in d.get('extra', []):
        st.markdown(f"##### {sec}")
        forms.section(sec)


#######################################################################################################################
# Inputs
#######################################################################################################################
def inputs(m):
    cfg = forms.cfg()
    with st.expander("Input files", expanded=False):
        files_form(cfg, upload=False)
    if m in vi.INPUTS:
        vi.INPUTS[m](cfg)
        return
    fields = active_files(cfg)
    if not fields:
        st.info("This mode reads no input workbooks.")
        return
    title = next((g['title'] for g, _, a in file_groups(cfg) if a and g['title']), None)
    if title:
        st.caption(f"Input files of the **{title}** (selected by the topology).")
    for tab, (var, kind, label, _) in zip(st.tabs([f"{f[2]} · {cfg[A.RUN_SECTION].get(f[0])}" for f in fields]),
                                          fields):
        with tab:
            vi.workbook(kind[5:], cfg[A.RUN_SECTION].get(var), key=f"m{m}|{var}")


#######################################################################################################################
# Run
#######################################################################################################################
def _summary(m):
    cfg = forms.cfg()
    G = cfg[A.RUN_SECTION]
    rows = [('Simulation', f"{G.get('Name', '')}  →  {os.path.relpath(A.results_root(cfg), A.ROOT)}/")]
    rows += [(label, str(G.get(var))) for var, kind, label, _ in active_files(cfg)]
    if A.MODES[m].get('sub'):
        rows.append(('Sub-mode', forms.option_label(*A.MODES[m]['sub'])))
    return rows


def run_view(m):
    cfg = forms.cfg()
    c1, c2 = st.columns([3, 2])
    with c1:
        with st.container(border=True):
            st.markdown("**Run settings**")
            for k, v in _summary(m):
                st.markdown(f"<span style='color:#52514e'>{k}</span> &nbsp; {v}", unsafe_allow_html=True)
    with c2:
        start = st.button(f"Run {A.MODES[m]['title']}", type='primary', icon=':material/play_arrow:',
                          width='stretch', key=f"run|{m}")
        st.caption("The simulation runs locally in this process; the console output is streamed below.")

    if start:
        from gui.core.runner import run_mode
        with st.status(f"Running {A.MODES[m]['title']} …", expanded=True) as status:
            box = st.empty()
            run = run_mode(cfg, m, forms.schema(), stream=lambda txt: box.code(txt[-8000:], language=None,
                                                                                height=380))
            if run['error']:
                status.update(label=f"Stopped with an error after {run['duration']:.1f} s", state='error')
            else:
                status.update(label=f"Finished in {run['duration']:.1f} s", state='complete', expanded=False)
        st.session_state.runs[m] = run
        try:
            files.save_run(run)
        except Exception as e:
            st.warning(f"The run could not be stored for the results browser: {e}")
        if not run['error']:
            st.session_state[f"goto|{m}"] = 'Results'
            st.rerun()
        else:
            st.error(run['error'])
            with st.expander("Console log"):
                st.code(run['log'][-20000:], language=None)


#######################################################################################################################
# Results
#######################################################################################################################
def results_view(m):
    cfg = forms.cfg()
    root = A.results_root(cfg)
    past = []
    if os.path.isdir(root):
        for f in files.list_runs('', root=root):
            meta = files.run_meta(f)
            if meta and meta.get('mode') == m:
                past.append((f, meta))
    opts = (['session'] if m in st.session_state.runs else []) + [f for f, _ in past]
    if not opts:
        vr.run_view(None)
        return
    lab = {f: f"{os.path.basename(f)} · {meta.get('time', '')}" for f, meta in past}
    lab['session'] = "Latest run of this session"
    sel = st.selectbox("Results", opts, format_func=lambda o: lab[o], key=f"res|pick|{m}",
                       help=f"GUI runs of this mode in {os.path.relpath(root, A.ROOT)}/")
    run = st.session_state.runs[m] if sel == 'session' else files.load_run(sel)
    vr.run_view(run, key=f"res|{m}")


#######################################################################################################################
# Page
#######################################################################################################################
def make_page(m):
    def page():
        d = A.MODES[m]
        if A.MODE_VAR:                                                  # activity rules and saved configs see the mode
            forms.cfg()[A.RUN_SECTION][A.MODE_VAR] = m
        st.title(f"{m} · {d['title']}" if len(A.MODES) > 1 else d['title'])
        st.caption(d['text'])
        vs = views(m)
        key = f"view|{m}"
        if f"goto|{m}" in st.session_state:
            st.session_state[key] = st.session_state.pop(f"goto|{m}")
        if st.session_state.get(key) not in vs:
            st.session_state[key] = vs[0]
        view = st.segmented_control("View", vs, key=key, label_visibility='collapsed',
                                    format_func=lambda v: f"{VIEW_ICONS[v]} {v}") or vs[0]
        st.divider()
        {'Workflow': workflow.view, 'Parameters': parameters, 'Inputs': inputs, 'Theory': theory_view,
         'Run': run_view, 'Results': results_view}[view](m)

    page.__name__ = f"mode{m}"
    return page
