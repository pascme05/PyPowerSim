#######################################################################################################################
#######################################################################################################################
# Title:        Power Electronics Simulation Toolkit (PyPowerSim)
# Topic:        Power Electronics
# File:         project.py
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
Setup page: run settings (formerly start.py, adapter.RUN_FIELDS), input files and the parameters of the run section
shared by all modes.
"""

#######################################################################################################################
# Import libs
#######################################################################################################################
from gui import adapter as A
from gui.core import files, forms

import os

import streamlit as st


#######################################################################################################################
# Widgets of the run settings
#######################################################################################################################
def _file(cfg, var, kind, label, help=None):
    opts = files.list_files(kind)
    cur = cfg[A.RUN_SECTION].get(var)
    if cur not in opts:
        opts = [cur] + opts
    cfg[A.RUN_SECTION][var] = st.selectbox(label, opts, index=opts.index(cur), key=f"g|{var}",
                                           help=help or f"Workbooks in the folder '{A.FOLDERS[kind]}/'.")


def run_field(cfg, var, kind, label, help=None):
    """Widget of one run setting, the value is written into the run section."""
    G = cfg[A.RUN_SECTION]
    if kind.startswith('expert:'):
        if not forms.expert():
            return
        kind = kind[7:]
    if kind == 'text':
        G[var] = st.text_input(label, value=str(G.get(var, '')), key=f"g|{var}", help=help)
    elif kind == 'toggle':
        v = st.toggle(label, value=bool(G.get(var)), key=f"g|{var}", help=help)
        G[var] = type(G.get(var, False))(v) if isinstance(G.get(var), (bool, int)) else v
    elif kind == 'int':
        G[var] = int(st.number_input(label, value=int(G.get(var, 0) or 0), step=1, key=f"g|{var}", help=help))
    elif kind.startswith('file:'):
        _file(cfg, var, kind[5:], label, help)
    elif kind.startswith('select:'):
        opts = kind[7:].split('|')
        cast = type(G.get(var)) if isinstance(G.get(var), (int, float)) else str
        opts = [cast(o) for o in opts]
        cur = G.get(var, opts[0])
        opts = opts if cur in opts else [cur] + opts
        G[var] = st.selectbox(label, opts, index=opts.index(cur), key=f"g|{var}", help=help)


def file_groups(cfg):
    """
    Input files grouped as in adapter.FILE_GROUPS: [(group dict, [RUN_FIELDS entries], active)]; one untitled group
    with all file settings if the adapter defines no groups.
    """
    fields = {f[0]: f for f in A.RUN_FIELDS if f[1].startswith('file:')}
    groups = getattr(A, 'FILE_GROUPS', None)
    if not groups:
        return [({'title': None}, list(fields.values()), True)] if fields else []
    out = []
    for g in groups:
        try:
            active = bool(g['active'](cfg))
        except Exception:
            active = True
        out.append((g, [fields[v] for v in g['fields'] if v in fields], active))
    return out


def active_files(cfg):
    """File settings used with the current configuration (all of them without FILE_GROUPS)."""
    return [f for _, fs, active in file_groups(cfg) if active for f in fs]


def _upload(cfg, fields, key):
    """Upload of a workbook into the folder of one file setting, which then uses it."""
    c1, c2 = st.columns([2, 3])
    with c1:
        var = st.selectbox("For", [f[0] for f in fields], key=f"up|for|{key}",
                           format_func=lambda v: next(f[2] for f in fields if f[0] == v))
    kind = next(f[1][5:] for f in fields if f[0] == var)
    with c2:
        up = st.file_uploader("Workbook", type=[e.lstrip('.') for e in A.FILE_EXT], key=f"up|{key}",
                              help=f"Stored in {A.FOLDERS[kind]}/ and selected for this setting.")
    if up is not None and st.session_state.get(f"upd|{key}") != up.file_id:
        st.session_state[f"upd|{key}"] = up.file_id
        cfg[A.RUN_SECTION][var] = files.save_upload(kind, up)
        st.session_state.pop(f"g|{var}", None)
        st.toast(f"{up.name} stored in {A.FOLDERS[kind]}/", icon=':material/upload:')
        st.rerun()


def _selector():
    """Parameter that selects the active file group (adapter.GROUP_SELECTOR), or None."""
    sel = getattr(A, 'GROUP_SELECTOR', None)
    try:
        return (sel[0], forms.find(*sel)) if sel else None
    except KeyError:
        return None


def files_form(cfg, upload=True):
    """Selection (and upload) of the input files named in RUN_FIELDS ('file:<kind>'), per FILE_GROUPS."""
    groups = file_groups(cfg)
    if not groups:
        return
    sel = _selector()
    if sel and groups[0][0]['title']:
        c1, c2 = st.columns([1, 3], vertical_alignment='bottom')
        with c1:
            forms.param(sel[0], sel[1], cfg)
        with c2:
            st.caption("The selected topology decides which group of input files the simulation uses; the other "
                       "group is kept for later.")
        groups = file_groups(cfg)                                     # the selector may have changed the group
    cols = st.columns(len(groups)) if len(groups) > 1 else [st.container()]
    for gi, ((g, fields, active), col) in enumerate(zip(groups, cols)):
        with col:
            with st.container(border=bool(g['title'])):
                if g['title']:
                    badge = (":green-badge[:material/check_circle: used by the current topology]" if active else
                             ":gray-badge[:material/pause_circle: not used by the current topology]")
                    st.markdown(f"##### {g.get('icon', '')} {g['title']} &nbsp; {badge}")
                    st.caption(g.get('text', ''))
                inner = st.columns(2 if len(fields) > 2 else max(len(fields), 1))
                for i, (var, kind, label, hlp) in enumerate(fields):
                    with inner[i % len(inner)]:
                        _file(cfg, var, kind[5:], label, hlp)
                if upload and fields:
                    with st.expander("Add a workbook", icon=':material/upload:'):
                        _upload(cfg, fields, key=str(gi))


#######################################################################################################################
# Page
#######################################################################################################################
def page():
    cfg = forms.cfg()
    if A.LOGO and os.path.isfile(os.path.join(A.ROOT, A.LOGO)):
        st.image(os.path.join(A.ROOT, A.LOGO), width=420)
    else:
        st.title(A.NAME)
    st.markdown(f"{A.TITLE} — {A.TAGLINE}")

    with st.container(border=True):
        st.markdown("#### :material/badge: Project")
        other = [f for f in A.RUN_FIELDS if not f[1].startswith('file:')]
        cols = st.columns(4)
        for i, (var, kind, label, hlp) in enumerate(other):
            with cols[i % 4]:
                run_field(cfg, var, kind, label, hlp)

    if any(f[1].startswith('file:') for f in A.RUN_FIELDS):
        with st.container(border=True):
            st.markdown("#### :material/description: Input files")
            files_form(cfg)

    if any(p['var'] not in A.RUN for p in forms.schema().get(A.RUN_SECTION, [])):
        with st.container(border=True):
            c1, c2 = st.columns([5, 1])
            with c1:
                st.markdown(f"#### :material/functions: {A.RUN_SECTION} parameters")
            with c2:
                forms.reset_section(A.RUN_SECTION)
            forms.section(A.RUN_SECTION, cols=4)

    st.markdown("#### Operating modes")
    cols = st.columns(min(len(A.MODES), 3))
    for i, (m, d) in enumerate(A.MODES.items()):
        with cols[i % len(cols)]:
            with st.container(border=True):
                st.markdown(f"**{d['icon']} {d['title']}**")
                st.caption(d['text'])
