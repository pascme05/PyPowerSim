#######################################################################################################################
#######################################################################################################################
# Title:        Power Electronics Simulation Toolkit (PyPowerSim)
# Topic:        Power Electronics
# File:         theory.py
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
Theory of the operating modes with interactive examples (view 'Theory' of the mode pages). The explorers start from
the current configuration and use the functions of the toolkit where possible: the PWM and the ideal time-domain model
of the topology classes (calcRef, calcPWM, calcTime), the loss interpolators of the device loader (Vce_2d, Eon_2d, …),
the thermal RC model (src/therm/calcTherRC.py), the THD function (src/general/helpFnc.py), the hysteresis controller
(src/cont/conHys.py) and the DAB power equation of src/general/genLoadInput.py. Changes in the explorers do not modify
the configuration.
"""

#######################################################################################################################
# Import libs
#######################################################################################################################
from gui.core import devices as D
from gui.core import toolkit as TK
from gui.core.style import show, line, color, fmt, table

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import streamlit as st

N_SWITCH = {'B2': 2, 'B4': 4, 'B6': 6, 'DAB': 8}


#######################################################################################################################
# Helper
#######################################################################################################################
def _intro(text):
    st.caption(text + " Values start from the current configuration; changes here do not modify it.")


def _dab(cfg):
    return cfg['Top']['sourceType'] == 'DAB'


def _prep(cfg):
    with st.spinner("Pre-processing with the toolkit …"):
        return TK.prepare(cfg)


def _spectrum(x):
    """Single-sided amplitude spectrum of one fundamental period (x sampled with endpoint): order -> peak value."""
    x = np.asarray(x, dtype=float)[:-1]
    X = np.abs(np.fft.rfft(x)) / len(x) * 2
    X[0] /= 2
    return np.arange(len(X)), X


def _thd(x, t):
    """THD ratio X_h,rms / X_1,rms of one period with the function of the toolkit (src/general/helpFnc.py)."""
    from src.general.helpFnc import thd
    return float(thd(np.asarray(x, dtype=float)[:-1], np.asarray(t, dtype=float)[:-1], 2 / np.sqrt(2), 1))


def _opts(sec, cat, var):
    from gui.core import forms
    try:
        return forms.find(sec, cat, var)['options'] or {}
    except KeyError:
        return {}


def _foster(para):
    """Foster network junction -> ambient of a switch (Rth_JC + Rth_CA) as arrays R [K/W], C [Ws/K]."""
    vec = para['Ther']['vec']
    R = np.concatenate([np.asarray(vec['Rth_JC'], float), np.asarray(vec['Rth_CA'], float)])
    C = np.concatenate([np.asarray(vec['Cth_JC'], float), np.asarray(vec['Cth_CA'], float)])
    return R, C, len(vec['Rth_JC'])


#######################################################################################################################
# Mode 0: modulation (inverter) and DAB power transfer
#######################################################################################################################
def _pwm(cfg, key):
    st.markdown("##### Pulse width modulation of a bridge leg")
    st.latex(r"v_{a0}(t) = s_a(t)\,\frac{V_{dc}}{2},\qquad \hat V_{a0,1} = M_i\,\frac{V_{dc}}{2}"
             r"\quad (M_i \le 1 \text{ sinusoidal},\ M_i \le 2/\sqrt{3} \text{ with zero sequence})")
    st.markdown("The switching function s_a = ±1 follows from the comparison of the modulation function x_a "
                "(reference + zero sequence n_0) with the carrier c. Its fundamental is proportional to the modulation "
                "index, the harmonics appear as side bands around multiples of the pulse number q = f_s / f_el. The "
                "load (R, L, back EMF) filters them, so the current is almost sinusoidal.")
    P = cfg['Par']['PWM']
    types = {k: v for k, v in _opts('Par', 'PWM', 'type').items() if k in ('FF', 'CB', 'SV')} or \
        {'FF': 'fundamental frequency', 'CB': 'carrier based', 'SV': 'space vector'}
    zeros = _opts('Par', 'PWM', 'zero') or {P['zero']: P['zero']}
    c = st.columns(4)
    with c[0]:
        Mi = st.slider("Modulation index Mi [–]", 0.05, 1.25, float(min(max(cfg['Dat']['stat']['Mi'], 0.05), 1.25)),
                       0.05, key=f"{key}|mi", help="Mi = 1: amplitude of the fundamental leg voltage = Vdc/2")
    with c[1]:
        typ = st.selectbox("PWM", list(types), index=list(types).index(P['type']) if P['type'] in types else 0,
                           format_func=lambda k: f"{k}: {types[k]}", key=f"{key}|typ")
    with c[2]:
        zero = st.selectbox("Zero sequence (CB)", list(zeros), index=list(zeros).index(P['zero'])
                            if P['zero'] in zeros else 0, format_func=lambda k: f"{k}: {zeros[k]}", key=f"{key}|zero",
                            disabled=typ != 'CB', help="Used by the carrier-based PWM (calcSeqCB)")
    with c[3]:
        q = st.number_input("Pulse number q = fs / fel [–]", 3, 99, int(min(max(round(P['fs'] / cfg['Top']['fel']), 3),
                                                                            99)),
                            key=f"{key}|q", help="Switching frequency fs = q·fel")
    fel = cfg['Top']['fel']
    fsim = float(min(cfg['Exp']['fsim'], 2e5))
    try:
        with st.spinner("calcRef · calcPWM · calcTime …"):
            d = TK.period(cfg, Mi, **{'Par.PWM.type': typ, 'Par.PWM.zero': zero, 'Par.PWM.fs': q * fel,
                                      'Exp.fsim': fsim})
    except Exception as e:
        st.warning(f"The toolkit could not compute this case ({type(e).__name__}: {e}).")
        return
    t = d['t'] * 1e3
    fig = make_subplots(rows=3, cols=1, shared_xaxes=True, vertical_spacing=0.06,
                        subplot_titles=("modulation function x_a, zero sequence n_0 and carrier c [p.u.]",
                                        "leg voltage v_a0 [V]", "phase current i_a [A]"))
    fig.add_trace(line(t, d['x']['A'], "x_a [p.u.]", 0), 1, 1)
    fig.add_trace(line(t, d['n0'], "n_0 [p.u.]", 2, width=1.5), 1, 1)
    fig.add_trace(line(t, d['c'], "carrier c [p.u.]", 7, width=1), 1, 1)
    fig.add_trace(line(t, d['ac']['v_a0'], "v_a0 [V]", 1, width=1), 2, 1)
    fig.add_trace(line(t, d['ac']['i_a'], "i_a [A]", 0), 3, 1)
    fig.update_xaxes(title_text="time t [ms]", row=3, col=1)
    show(fig, height=620, key=f"{key}|pwm", legend_below=True)
    h, Va = _spectrum(d['ac']['v_a0'])
    V1, Vdc = Va[1], d['Vdc']
    k = st.columns(4)
    k[0].metric("Fundamental v_a0 (peak)", fmt(V1, 'V'), help="Amplitude of the fundamental (FFT of one period)")
    k[1].metric("Linear theory Mi · Vdc/2", fmt(Mi * Vdc / 2, 'V'), help="Valid in the linear range of the PWM")
    k[2].metric("THD i_a", fmt(100 * _thd(d['ac']['i_a'], d['t']), '%'),
                help="RMS of the harmonics / RMS of the fundamental (src/general/helpFnc.thd)")
    k[3].metric("THD v_a0", fmt(100 * _thd(d['ac']['v_a0'], d['t']), '%'))
    m = (h >= 1) & (h <= 4 * q)
    fig = go.Figure(go.Bar(x=h[m], y=Va[m], marker_color=color(1), name="v_a0", width=0.7,
                           hovertemplate="h = %{x}<br>%{y:.4g} V<extra></extra>"))
    fig.update_layout(xaxis_title="harmonic order h = f / fel [–]", yaxis_title="amplitude of v_a0 [V]",
                      showlegend=False)
    fig.update_yaxes(type='log')
    show(fig, height=300, key=f"{key}|spec")
    big = np.argsort(Va[m])[::-1][:8]
    table(pd.DataFrame({'order h [–]': h[m][big], 'frequency [Hz]': h[m][big] * fel, 'amplitude [V]': Va[m][big],
                        'of fundamental [%]': 100 * Va[m][big] / V1 if V1 else np.nan}), hide_index=True)
    st.caption(f"Ideal switches, one period at fsim = {fmt(fsim, 'Hz')} (reduced for interactivity). Side bands "
               f"around h = q = {q}, 2q, 3q …; with a zero sequence (SV, THIPWM, DPWM) the triplen harmonics cancel "
               f"in the line-to-line voltages.")


def _mi_curve(cfg, key):
    st.markdown("##### Fundamental over the modulation index")
    st.markdown("In the linear range the fundamental of the leg voltage grows with Mi. Above the limit of the PWM "
                "method the pulses saturate (overmodulation) until the square wave (six-step) limit 4/π is reached. "
                "The points are computed with the toolkit (calcPWM, calcTime) for the PWM of the configuration.")
    P = cfg['Par']['PWM']
    Mis = np.round(np.arange(0.1, 1.31, 0.1), 2)
    pts = []
    try:
        with st.spinner(f"{len(Mis)} periods with the toolkit …"):
            for Mi in Mis:
                d = TK.period(cfg, Mi, **{'Exp.fsim': 1e5})
                pts.append(_spectrum(d['ac']['v_a0'])[1][1] / (d['Vdc'] / 2))
    except Exception as e:
        st.warning(f"The toolkit could not compute the curve ({type(e).__name__}: {e}).")
        return
    fig = go.Figure()
    x = np.linspace(0, 1.3, 50)
    fig.add_trace(line(x, x, "linear theory V̂1 = Mi·Vdc/2", 7, dash='dash'))
    fig.add_trace(line([0, 1.3], [4 / np.pi] * 2, "six-step limit 4/π", 3, dash='dot'))
    fig.add_trace(line(Mis, pts, f"toolkit: {P['type']}" + (f" / {P['zero']}" if P['type'] == 'CB' else ''), 0,
                       mode='lines+markers'))
    fig.add_vline(x=1, line_dash='dot', line_color='#888')
    fig.add_vline(x=2 / np.sqrt(3), line_dash='dot', line_color='#888')
    fig.update_layout(xaxis_title="modulation index Mi [–]", yaxis_title="fundamental V̂_a0,1 / (Vdc/2) [–]",
                      hovermode='x unified')
    show(fig, height=380, key=f"{key}|mi", legend_below=True)
    table(pd.DataFrame({'Mi [–]': Mis, 'V̂1 / (Vdc/2) [–]': pts, 'deviation from linear [%]':
                        100 * (np.array(pts) - Mis) / Mis}), hide_index=True)
    st.caption("Vertical lines: Mi = 1 (limit of sinusoidal PWM) and Mi = 2/√3 ≈ 1.155 (limit with zero sequence / "
               "space vector PWM).")


def _dab_power(cfg, key):
    st.markdown("##### Power transfer of the dual active bridge (single phase shift)")
    st.latex(r"P = \frac{N\,V_1\,V_2}{\omega L_k}\left(\varphi - \frac{\varphi^2}{\pi}\right),\qquad "
             r"\omega = 2\pi f_s,\quad L_k = L_p + N^2 L_s,\quad N = N_p/N_s")
    st.markdown("Both bridges generate square-wave voltages; the phase shift φ between them drives the current through "
                "the leakage inductance L_k. The power is maximal at φ = π/2. Zero-voltage switching (ZVS) requires "
                "the current to have the right sign at the switching instants, which depends on the voltage ratio "
                "d = N·V2/V1 and on φ (power equation of src/general/genLoadInput.py).")
    try:
        s, para, _, _ = _prep(cfg)
        N0, Lk0 = para['Tra']['Elec']['con']['N'], para['Tra']['Elec']['con']['Lk']
    except Exception:
        s, N0, Lk0 = cfg, 1.0, 50e-6
    c = st.columns(5)
    with c[0]:
        V1 = st.number_input("V1 input [V]", 1.0, 5000.0, float(abs(s['Dat']['stat']['Vdc'])), key=f"{key}|v1")
    with c[1]:
        V2 = st.number_input("V2 output [V]", 1.0, 5000.0, float(abs(s['Dat']['stat'].get('Vo', V1 / N0))),
                             key=f"{key}|v2")
    with c[2]:
        N = st.number_input("N = Np/Ns [–]", 0.05, 50.0, float(N0), key=f"{key}|n")
    with c[3]:
        fs = st.number_input("fs [kHz]", 1.0, 2000.0, float(cfg['Par']['PWM']['fs']) / 1e3, key=f"{key}|fs") * 1e3
    with c[4]:
        Lk = st.number_input("Lk [µH]", 0.1, 10000.0, float(Lk0) * 1e6, key=f"{key}|lk") * 1e-6
    phi0 = float(np.degrees(abs(s['Dat']['stat'].get('PhiDAB', np.radians(30)))))
    phi_d = st.slider("Phase shift φ [°]", 0.0, 90.0, float(min(phi0, 90.0)), 1.0, key=f"{key}|phi")
    phi = np.radians(phi_d)
    w = 2 * np.pi * fs
    V2r = N * V2                                                      # output voltage referred to the primary
    p = lambda f: N * V1 * V2 / (w * Lk) * (f - f ** 2 / np.pi)
    th = np.linspace(0, np.pi, 400)
    i0 = -(V1 * np.pi + V2r * (2 * phi - np.pi)) / (2 * w * Lk)
    i_phi = (V1 * (2 * phi - np.pi) + V2r * np.pi) / (2 * w * Lk)
    i = np.where(th < phi, i0 + (V1 + V2r) * th / (w * Lk), i_phi + (V1 - V2r) * (th - phi) / (w * Lk))
    I_rms = float(np.sqrt(np.mean(i ** 2)))
    k = st.columns(5)
    k[0].metric("Power P", fmt(p(phi), 'W'))
    k[1].metric("Maximum (φ = 90°)", fmt(p(np.pi / 2), 'W'))
    k[2].metric("Inductor current (RMS, prim.)", fmt(I_rms, 'A'))
    k[3].metric("ZVS primary bridge", "yes" if i0 < 0 else "no", help="i(0) < 0 when the primary bridge switches")
    k[4].metric("ZVS secondary bridge", "yes" if i_phi > 0 else "no", help="i(φ) > 0 when the secondary switches")
    c1, c2 = st.columns(2)
    with c1:
        f = np.linspace(0, np.pi / 2, 100)
        fig = go.Figure([line(np.degrees(f), p(f) / 1e3, "P(φ) [kW]", 0),
                         go.Scatter(x=[phi_d], y=[p(phi) / 1e3], mode='markers', marker=dict(size=11, color=color(1)),
                                    name="operating point")])
        fig.update_layout(xaxis_title="phase shift φ [°]", yaxis_title="power P [kW]")
        show(fig, height=340, key=f"{key}|pphi", legend_below=True)
    with c2:
        ang = np.concatenate([th, th + np.pi])
        tt = ang / w * 1e6
        ii = np.concatenate([i, -i])
        v1 = np.where(ang < np.pi, V1, -V1)
        v2 = np.where((ang >= phi) & (ang < phi + np.pi), V2r, -V2r)
        fig = make_subplots(specs=[[{'secondary_y': True}]])
        fig.add_trace(line(tt, ii, "i_L (primary) [A]", 0), secondary_y=False)
        fig.add_trace(line(tt, v1, "v_1 [V]", 1, width=1), secondary_y=True)
        fig.add_trace(line(tt, v2, "N·v_2 [V]", 2, width=1, dash='dot'), secondary_y=True)
        fig.update_layout(xaxis_title="time t [µs]")
        fig.update_yaxes(title_text="current [A]", secondary_y=False)
        fig.update_yaxes(title_text="voltage [V]", secondary_y=True, showgrid=False)
        show(fig, height=340, key=f"{key}|wave", legend_below=True)
    st.caption("Ideal, lossless model with piecewise linear current (one switching period). The toolkit simulates the "
               "same circuit in the time domain incl. dead time, device and transformer losses.")


def mode0(cfg):
    _intro("The sweep varies the operating point (modulation index or phase shift) and compares the simulated RMS "
           "values and distortions with analytical solutions.")
    if _dab(cfg):
        t1, = st.tabs(["DAB power transfer"])
        with t1:
            _dab_power(cfg, 'th0|dab')
        return
    t1, t2 = st.tabs(["PWM explorer", "Fundamental over Mi"])
    with t1:
        _pwm(cfg, 'th0|pwm')
    with t2:
        _mi_curve(cfg, 'th0|mi')


#######################################################################################################################
# Mode 1: losses and thermal ripple
#######################################################################################################################
def _device(cfg):
    """(para of the switch, setup, name) of the switches of the configuration (primary bridge for the DAB)."""
    s, para, _, _ = _prep(cfg)
    key = 'SwiPri' if _dab(cfg) else 'Swi'
    return para[key], s, cfg['Exp'][key]


def _interp(para, name, Tj, x):
    """Interpolator of the toolkit (para['Elec']['tab'][name + '_2d']) or the constant value."""
    f = para['Elec']['tab'].get(f"{name}_2d")
    x = np.abs(np.asarray(x, dtype=float))
    if f is not None:
        return np.nan_to_num(np.asarray(f((np.full_like(x, Tj), x)), dtype=float))
    return np.full_like(x, float(para['Elec']['con'].get(name, 0.0) or 0.0))


def _losses(cfg, key):
    st.markdown("##### Conduction and switching losses of one switch over a fundamental period")
    st.latex(r"P_{T,c} = \frac{1}{2\pi}\int_{i>0} d(\theta)\,v_{CE}(i)\,i\,d\theta,\qquad "
             r"P_{T,s} = \frac{f_s}{2\pi}\int_{i>0}\big(E_{on}(i)+E_{off}(i)\big)\frac{V_{dc}}{V_{nom}}\,d\theta")
    st.latex(r"P_{D,c} = \frac{1}{2\pi}\int_{i<0} d(\theta)\,v_{F}(|i|)\,|i|\,d\theta,\qquad "
             r"P_{D,s} = \frac{f_s}{2\pi}\int_{i<0}E_{rec}(|i|)\frac{V_{dc}}{V_{nom}}\,d\theta,\qquad "
             r"d(\theta) = \tfrac12\big(1 + M_i\sin\theta\big)")
    st.markdown("The explorer averages the loss models of the device workbook (the interpolators Vce_2d, Vfd_2d, "
                "Eon_2d, Eoff_2d, Erec_2d of the toolkit loader) over a sinusoidal current i = Î·sin(θ − φ) with the duty "
                "cycle of sinusoidal PWM; while the upper switch is turned on, its transistor carries i > 0 and its "
                "anti-parallel diode i < 0. The simulation evaluates the same models per sample and switching event, so "
                "its results differ slightly (zero sequence, dead time, ripple current).")
    try:
        para, s, name = _device(cfg)
    except Exception as e:
        st.warning(f"The device could not be loaded ({type(e).__name__}: {e}).")
        return None
    st.caption(f"Device: **{name}** ({s['Par']['Elec']['SwiType']}), Vnom = "
               f"{fmt(para['Elec']['con'].get('Vnom'), 'V')} (reference voltage of the switching energies)")
    I0 = float(abs(s['Dat']['stat'].get('Io', 20)) * np.sqrt(2))
    phi0 = float(np.degrees(abs(s['Dat']['stat'].get('PhiVI', 0.0))))
    c = st.columns(6)
    with c[0]:
        I = st.number_input("Î peak current [A]", 0.1, 5000.0, round(max(I0, 0.1), 1), key=f"{key}|i")
    with c[1]:
        Vdc = st.number_input("Vdc [V]", 1.0, 5000.0, float(abs(s['Dat']['stat']['Vdc'])), key=f"{key}|vdc")
    with c[2]:
        fs = st.number_input("fs [kHz]", 0.01, 2000.0, float(s['Par']['PWM']['fs']) / 1e3, key=f"{key}|fs") * 1e3
    with c[3]:
        Mi = st.number_input("Mi [–]", 0.0, 1.2, float(min(abs(s['Dat']['stat']['Mi']), 1.2)), 0.05, key=f"{key}|mi")
    with c[4]:
        phi = np.radians(st.number_input("φ load angle [°]", -90.0, 90.0, round(phi0, 1), key=f"{key}|phi",
                                         help="angle between voltage and current (PhiVI of the operating point)"))
    with c[5]:
        Tj = st.number_input("Tj [°C]", -40.0, 200.0, float(s['Dat']['stat'].get('Tj', 25)), key=f"{key}|tj")
    th = np.linspace(0, 2 * np.pi, 1441)[:-1]
    i = I * np.sin(th - phi)
    d = 0.5 * (1 + Mi * np.sin(th))
    Vnom = float(para['Elec']['con'].get('Vnom') or Vdc)
    pos, neg = i > 0, i < 0
    vT, vD = _interp(para, 'Vce', Tj, i), _interp(para, 'Vfd', Tj, i)
    E_on, E_off, E_rec = (_interp(para, n, Tj, i) for n in ('Eon', 'Eoff', 'Erec'))
    res = {'transistor conduction': np.mean(np.where(pos, d * vT * np.abs(i), 0)),
           'transistor switching': fs * np.mean(np.where(pos, (E_on + E_off) * Vdc / Vnom, 0)),
           'diode conduction': np.mean(np.where(neg, d * vD * np.abs(i), 0)),     # gate on, current reversed
           'diode switching': fs * np.mean(np.where(neg, E_rec * Vdc / Vnom, 0))}
    total = sum(res.values())
    n = N_SWITCH.get(cfg['Top']['sourceType'], 6)
    k = st.columns(4)
    k[0].metric("Losses per switch", fmt(total, 'W'), help="transistor + diode of one switch position")
    k[1].metric(f"Bridge ({n} switches)", fmt(n * total, 'W'))
    k[2].metric("Switching share", fmt(100 * (res['transistor switching'] + res['diode switching']) / total, '%')
                if total else '–')
    k[3].metric("Transistor / diode", f"{fmt(res['transistor conduction'] + res['transistor switching'], 'W')} / "
                                      f"{fmt(res['diode conduction'] + res['diode switching'], 'W')}")
    c1, c2 = st.columns(2)
    with c1:
        fig = go.Figure(go.Bar(x=list(res), y=list(res.values()), marker_color=[color(j) for j in range(4)],
                               text=[fmt(v, 'W') for v in res.values()], textposition='outside',
                               hovertemplate="%{x}<br>%{y:.4g} W<extra></extra>"))
        fig.update_layout(yaxis_title="mean loss P [W]", showlegend=False)
        show(fig, height=340, key=f"{key}|bar")
    with c2:
        x = np.linspace(0, max(I, 1) * 1.2, 120)
        fig = make_subplots(specs=[[{'secondary_y': True}]])
        fig.add_trace(line(x, _interp(para, 'Vce', Tj, x), "v_CE transistor [V]", 0), secondary_y=False)
        fig.add_trace(line(x, _interp(para, 'Vfd', Tj, x), "v_F diode [V]", 1), secondary_y=False)
        fig.add_trace(line(x, 1e3 * (_interp(para, 'Eon', Tj, x) + _interp(para, 'Eoff', Tj, x)) * Vdc / Vnom,
                           "E_on + E_off at Vdc [mJ]", 2, dash='dash'), secondary_y=True)
        fig.add_trace(line(x, 1e3 * _interp(para, 'Erec', Tj, x) * Vdc / Vnom, "E_rec at Vdc [mJ]", 3, dash='dash'),
                      secondary_y=True)
        fig.add_vrect(x0=0, x1=I, fillcolor='#2a78d6', opacity=0.06, line_width=0)
        fig.update_layout(xaxis_title="current |i| [A]", hovermode='x unified')
        fig.update_yaxes(title_text="on-state voltage [V]", secondary_y=False)
        fig.update_yaxes(title_text="switching energy [mJ]", secondary_y=True, showgrid=False)
        show(fig, height=340, key=f"{key}|curves", legend_below=True)
    st.caption(f"Shaded: current range of the operating point (0 … Î). Device models at Tj = {Tj:g} °C, energies "
               f"scaled linearly from Vnom to Vdc as in src/elec/calcLossSwi.py.")
    return total, res


def _ripple(cfg, key, P0=None):
    st.markdown("##### Junction temperature ripple over a fundamental period")
    st.latex(r"T_j(t) = T_a + \sum_i \Delta T_i(t),\qquad "
             r"\tau_i\,\frac{d\Delta T_i}{dt} + \Delta T_i = R_{th,i}\,P_v(t),\qquad \tau_i = R_{th,i}C_{th,i}")
    st.markdown("A switch conducts only during one half wave, so its losses pulsate with the fundamental frequency. "
                "Branches of the Foster network with time constants shorter than the period follow the losses and "
                "cause a temperature ripple ΔTj; the slow branches only see the mean loss. Low output frequencies "
                "(e.g. start-up of a motor) therefore give large temperature cycles. Computed with "
                "src/therm/calcTherRC.py.")
    from src.therm.calcTherRC import calcTherRC
    try:
        para, s, name = _device(cfg)
    except Exception as e:
        st.warning(f"The device could not be loaded ({type(e).__name__}: {e}).")
        return
    R, C, _ = _foster(para)
    c = st.columns(3)
    with c[0]:
        P = st.number_input("Mean loss of the transistor P [W]", 0.1, 5000.0, float(round(P0 or 30.0, 1)),
                            key=f"{key}|p")
    with c[1]:
        fel = st.number_input("Fundamental frequency fel [Hz]", 0.1, 1e5, float(cfg['Top']['fel']), key=f"{key}|fel")
    with c[2]:
        Ta = st.number_input("Ambient / coolant Ta [°C]", -40.0, 150.0, float(s['Dat']['trans'].get('Tc', 25)),
                             key=f"{key}|ta")

    def tj(f, periods=3, n=400):
        t = np.linspace(0, periods / f, periods * n + 1)
        th = 2 * np.pi * f * t
        pv = np.where(np.sin(th) > 0, 4 * P * np.sin(th) ** 2, 0.0)            # half wave, mean P
        T, _ = calcTherRC(R * P, Ta, pv, t, R, C)                             # start at the mean temperature rise
        return t, T

    t, T = tj(fel)
    last = T[-401:]
    k = st.columns(4)
    k[0].metric("Mean Tj", fmt(Ta + P * R.sum(), '°C'), help="Ta + P·ΣRth")
    k[1].metric("Ripple ΔTj", fmt(last.max() - last.min(), 'K'), help="max − min in the last period")
    k[2].metric("Peak Tj", fmt(last.max(), '°C'))
    k[3].metric("Rth junction → ambient", fmt(R.sum(), 'K/W'))
    c1, c2 = st.columns(2)
    with c1:
        fig = go.Figure([line(t * 1e3, T, "T_j [°C]", 1)])
        fig.update_layout(xaxis_title="time t [ms]", yaxis_title="junction temperature Tj [°C]")
        show(fig, height=320, key=f"{key}|tj")
    with c2:
        fs_ = np.logspace(-1, 3, 17)
        dT = []
        for f in fs_:
            _, T2 = tj(f, periods=2, n=200)
            dT.append(T2[-201:].max() - T2[-201:].min())
        fig = go.Figure([line(fs_, dT, "ΔTj [K]", 0, mode='lines+markers'),
                         go.Scatter(x=[fel], y=[last.max() - last.min()], mode='markers',
                                    marker=dict(size=11, color=color(1)), name="current fel")])
        fig.update_layout(xaxis_title="fundamental frequency fel [Hz]", yaxis_title="temperature ripple ΔTj [K]")
        fig.update_xaxes(type='log')
        show(fig, height=320, key=f"{key}|dtf", legend_below=True)
    st.caption(f"Device {name}: Foster network junction → case (Rth_JC, Cth_JC) and case → ambient (Rth_CA, Cth_CA) "
               f"of the workbook; loss pulse 4P·sin²(θ) during the positive half wave (mean P).")


def mode1(cfg):
    _intro("The steady state simulates cyc fundamental periods with the electrical, loss and thermal models and "
           "iterates until the temperatures change less than the tolerance.")
    if _dab(cfg):                                                     # loss explorers assume a sinusoidal current
        t1, t2 = st.tabs(["DAB power transfer", "Steady temperature of a switch"])
        with t1:
            _dab_power(cfg, 'th1|dab')
        with t2:
            st.caption("The switches of the DAB conduct square-wave currents at the switching frequency; the loss "
                       "pulses are much faster than the thermal time constants, so the temperature ripple is small "
                       "and the steady temperature is Ta + P·ΣRth (end of the curve).")
            _heating(cfg, 'th1|heat')
        return
    t1, t2 = st.tabs(["Losses of a switch", "Temperature ripple"])
    with t1:
        out = _losses(cfg, 'th1|loss')
    with t2:
        _ripple(cfg, 'th1|rip', P0=(out[1]['transistor conduction'] + out[1]['transistor switching']) if out else None)


#######################################################################################################################
# Mode 2: transients
#######################################################################################################################
def _heating(cfg, key):
    st.markdown("##### Heating of a device after a load step")
    st.latex(r"T_j(t) = T_a + P\,Z_{th,ja}(t),\qquad Z_{th,ja}(t) = \sum_i R_{th,i}\left(1-e^{-t/\tau_i}\right)")
    st.markdown("The transient simulation starts at the initial temperatures (Dat › trans: Tj, Tc) and follows the "
                "heating until tmax. The fast branches of the Foster network (chip, solder) settle within "
                "milliseconds, the case and heatsink branches need seconds to minutes — a transient run shorter "
                "than the largest time constant does not reach the steady state.")
    from src.therm.calcTherRC import calcTherRC
    try:
        para, s, name = _device(cfg)
    except Exception as e:
        st.warning(f"The device could not be loaded ({type(e).__name__}: {e}).")
        return
    R, C, n_jc = _foster(para)
    tmax = float(s['Dat']['trans']['tmax'])
    c = st.columns(3)
    with c[0]:
        P = st.number_input("Loss step P [W]", 0.1, 5000.0, 30.0, key=f"{key}|p")
    with c[1]:
        Ta = st.number_input("Initial / ambient temperature [°C]", -40.0, 150.0,
                             float(s['Dat']['trans'].get('Tj', 25)), key=f"{key}|ta")
    with c[2]:
        tm = st.number_input("tmax [s]", 1e-4, 1e5, tmax, format='%g', key=f"{key}|tm",
                             help="Simulated time of the transient mode (Dat › trans › tmax)")
    tau = R * C
    t = np.concatenate([[0], np.logspace(np.log10(max(tau.min() / 20, 1e-7)), np.log10(max(tau.max() * 10, tm)),
                                         600)])
    T, _ = calcTherRC(np.zeros(len(R)), Ta, np.full(len(t), P), t, R, C)
    reach = float(np.interp(tm, t, T))
    final = Ta + P * R.sum()
    k = st.columns(4)
    k[0].metric("Final temperature", fmt(final, '°C'), help="Ta + P·ΣRth (steady state)")
    k[1].metric(f"At tmax = {fmt(tm, 's')}", fmt(reach, '°C'))
    k[2].metric("Heating reached", fmt(100 * (reach - Ta) / (final - Ta), '%') if final > Ta else '–')
    k[3].metric("Largest time constant", fmt(tau.max(), 's'))
    fig = go.Figure([line(t[1:], T[1:], "T_j (calcTherRC) [°C]", 1)])
    fig.add_trace(line(t[1:], Ta + P * D.zth(R, C, t[1:]), "Ta + P·Zth(t) (analytical) [°C]", 0, dash='dot'))
    fig.add_vline(x=tm, line_dash='dash', line_color='#888', annotation_text="tmax")
    fig.update_layout(xaxis_title="time t [s]", yaxis_title="junction temperature Tj [°C]", hovermode='x unified')
    fig.update_xaxes(type='log')
    show(fig, height=360, key=f"{key}|heat", legend_below=True)
    table(pd.DataFrame({'branch': [f"{'JC' if j < n_jc else 'CA'} {j + 1}" for j in range(len(R))],
                        'Rth [K/W]': R, 'Cth [Ws/K]': C, 'τ = Rth·Cth [s]': tau,
                        'settled at tmax [%]': 100 * (1 - np.exp(-tm / tau))}), hide_index=True)
    st.caption(f"Device {name}: junction → case (JC) and case → ambient (CA) Foster networks of the workbook. "
               f"A branch has settled to 63 % after τ and to 99 % after 4.6 τ.")


def _startup(cfg, key):
    st.markdown("##### Electrical start-up of the RL load")
    st.latex(r"i(t) = \hat I\left[\sin(\omega t + \alpha - \varphi) - \sin(\alpha - \varphi)\,e^{-t/\tau}\right],"
             r"\qquad \hat I = \frac{\hat V}{\sqrt{R^2 + (\omega L)^2}},\quad \tau = L/R,\quad"
             r"\varphi = \arctan(\omega L / R)")
    st.markdown("When the fundamental voltage is switched on, the current contains a decaying DC component that "
                "depends on the switching angle α. It decays with τ = L/R; the transient analysis shows it in the first "
                "periods (the PWM ripple is added on top).")
    R, L, fel = float(cfg['Top']['R']), float(cfg['Top']['L']), float(cfg['Top']['fel'])
    Vdc, Mi = float(abs(cfg['Dat']['stat']['Vdc'])), float(abs(cfg['Dat']['stat']['Mi']))
    c = st.columns(4)
    with c[0]:
        R = st.number_input("R [Ω]", 0.001, 1e4, max(R, 0.001), key=f"{key}|r")
    with c[1]:
        L = st.number_input("L [mH]", 0.0, 1e4, L * 1e3, key=f"{key}|l") * 1e-3
    with c[2]:
        V = st.number_input("V̂ fundamental [V]", 0.0, 1e4, Mi * Vdc / 2, key=f"{key}|v", help="Mi·Vdc/2")
    with c[3]:
        a = np.radians(st.slider("Switching angle α [°]", 0, 180, 0, 5, key=f"{key}|a"))
    w = 2 * np.pi * fel
    Z, ph, tau = np.hypot(R, w * L), np.arctan2(w * L, R), L / R
    t = np.linspace(0, max(5 / fel, 5 * tau), 2000)
    I = V / Z
    dc = -I * np.sin(a - ph) * np.exp(-t / tau) if tau > 0 else np.zeros_like(t)
    i = I * np.sin(w * t + a - ph) + dc
    k = st.columns(3)
    k[0].metric("τ = L/R", fmt(tau, 's'))
    k[1].metric("Steady amplitude Î", fmt(I, 'A'))
    k[2].metric("Peak during start-up", fmt(np.max(np.abs(i)), 'A'))
    fig = go.Figure([line(t * 1e3, i, "i(t) [A]", 0), line(t * 1e3, dc, "DC component [A]", 1, dash='dash')])
    fig.update_layout(xaxis_title="time t [ms]", yaxis_title="current i [A]", hovermode='x unified')
    show(fig, height=320, key=f"{key}|i", legend_below=True)


def mode2(cfg):
    _intro("The transient analysis simulates tmax seconds with the coupled electrical, loss and thermal models, "
           "starting from the initial temperatures.")
    t1, t2 = st.tabs(["Heating of the devices", "Electrical start-up"])
    with t1:
        _heating(cfg, 'th2|heat')
    with t2:
        if _dab(cfg):
            st.info("The start-up explorer is shown for the inverter topologies (RL load).")
        else:
            _startup(cfg, 'th2|st')


#######################################################################################################################
# Mode 3: closed loop
#######################################################################################################################
def _hysteresis(cfg, key):
    st.markdown("##### Hysteresis current control")
    st.latex(r"s_a = \begin{cases} -1 & i - i_{ref} \ge \delta \\ +1 & i - i_{ref} < -\delta \\ s_a & \text{else}"
             r"\end{cases},\qquad \delta = \frac{h}{100}\,\max|i_{ref}|,\qquad "
             r"L\frac{di}{dt} = s_a\frac{V_{dc}}{2} - R\,i")
    st.markdown("The controller of the toolkit (src/cont/conHys.py) compares the current with the reference once per "
                "control period 1/fc and holds the switching state in between. Narrow bands track the reference "
                "closely but switch more often; the sampling with fc limits the switching frequency and lets the "
                "current overshoot the band. As in the toolkit the reference starts at 50 % and steps to 100 % at "
                "half of the simulated time. The explorer uses one bridge leg with the RL load (half-bridge "
                "equivalent); in the B6 bridge the phases are coupled through the star point.")
    from src.cont.conHys import conHys
    Cn = cfg['Par']['Cont']
    s = cfg
    try:
        s = _prep(cfg)[0]
    except Exception:
        pass
    I0 = float(abs(s['Dat']['stat'].get('Io', 20)) * np.sqrt(2))
    c = st.columns(5)
    with c[0]:
        h = st.number_input("Band h [%]", 0.1, 50.0, float(Cn['hys']), key=f"{key}|h")
    with c[1]:
        fc = st.number_input("fc [kHz]", 0.5, 500.0, float(Cn['fc']) / 1e3, key=f"{key}|fc") * 1e3
    with c[2]:
        Iref = st.number_input("Reference Î [A]", 0.1, 5000.0, round(I0, 1), key=f"{key}|i")
    with c[3]:
        L = st.number_input("L [mH]", 0.01, 1e4, max(float(cfg['Top']['L']) * 1e3, 0.01), key=f"{key}|l") * 1e-3
    with c[4]:
        Vdc = st.number_input("Vdc [V]", 1.0, 5000.0, float(abs(cfg['Dat']['stat']['Vdc'])), key=f"{key}|vdc")
    R, fel = max(float(cfg['Top']['R']), 1e-3), float(cfg['Top']['fel'])
    T = 2 / fel
    n_sub = 10                                                        # integration steps per control period
    n_con = int(np.ceil(T * fc))
    dt = 1 / fc / n_sub
    a = np.exp(-R / L * dt)
    i_v, r_v, s_v = np.zeros(n_con * n_sub), np.zeros(n_con * n_sub), np.zeros(n_con * n_sub)
    s_act, i_act = np.ones(n_sub), 0.0
    tol_base = Iref * h / 100
    for k in range(n_con):
        tk = k / fc
        scale = 0.5 if k / n_con < 0.5 else 1.0
        s_act, _ = conHys(i_act, scale * Iref * np.sin(2 * np.pi * fel * tk), s_act, tol_base * scale)
        for j in range(n_sub):
            i_act = a * i_act + (1 - a) * s_act[0] * Vdc / 2 / R      # exact solution of the RL load per step
            n = k * n_sub + j
            i_v[n], s_v[n] = i_act, s_act[0]
            r_v[n] = scale * Iref * np.sin(2 * np.pi * fel * (tk + (j + 1) * dt))
    t = (np.arange(len(i_v)) + 1) * dt
    sw = int(np.sum(np.abs(np.diff(s_v)) > 0))
    err = i_v - r_v
    half = len(t) // 2
    k = st.columns(4)
    k[0].metric("Average switching frequency", fmt(sw / 2 / T, 'Hz'), help="transitions / 2 / simulated time")
    k[1].metric("Band ±δ (100 % reference)", fmt(tol_base, 'A'))
    k[2].metric("Max |error| (2nd half)", fmt(np.max(np.abs(err[half:])), 'A'))
    k[3].metric("RMS error (2nd half)", fmt(np.sqrt(np.mean(err[half:] ** 2)), 'A'))
    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, row_heights=[0.75, 0.25], vertical_spacing=0.05)
    fig.add_trace(line(t * 1e3, i_v, "current i [A]", 0, width=1.2), 1, 1)
    fig.add_trace(line(t * 1e3, r_v, "reference i_ref [A]", 1, width=1.5), 1, 1)
    band_v = tol_base * np.where(np.arange(len(t)) < half, 0.5, 1.0)
    fig.add_trace(line(t * 1e3, r_v + band_v, "band i_ref ± δ [A]", 7, width=0.8, dash='dot'), 1, 1)
    fig.add_trace(line(t * 1e3, r_v - band_v, "band", 7, width=0.8, dash='dot', showlegend=False), 1, 1)
    fig.add_trace(line(t * 1e3, s_v, "switching state s_a [–]", 2, width=1, line_shape='hv'), 2, 1)
    fig.update_xaxes(title_text="time t [ms]", row=2, col=1)
    fig.update_yaxes(title_text="current [A]", row=1, col=1)
    fig.update_yaxes(title_text="s_a [–]", row=2, col=1)
    show(fig, height=520, key=f"{key}|hys", legend_below=True)
    st.caption("In B2, B4 and B6 the setting PI also uses the hysteresis control (calcCON of the topology classes); "
               "Kp and Ki act only in the DAB.")


def _dab_pi(cfg, key):
    st.markdown("##### PI control of the DAB output current")
    st.latex(r"\varphi[k] = \operatorname{sat}\!\left(K_p\,e[k] + K_i \sum_{j\le k} e[j]\,T_c\right),\qquad "
             r"e = i_{ref} - i_{act},\quad T_c = 1/f_c,\quad |\varphi| \le \pi/2")
    st.markdown("For the DAB the PI controller sets the phase shift once per control period (classDAB.calcCON): the "
                "output is limited by satMin/satMax and ±π/2, the integrator stops while the output is saturated "
                "(anti-windup). The static gain dP/dφ of the power stage follows from the power equation below — it "
                "decreases towards φ = 90°, so a fixed PI tuning reacts slower at high load. With the setting HY the "
                "DAB uses the hysteresis control of the inverters instead.")
    _dab_power(cfg, f"{key}|p")


def mode3(cfg):
    _intro("The closed loop simulates tmax·fc control periods; in each period the controller updates the switching "
           "state (inverter: hysteresis) or the phase shift (DAB: PI) from the measured current.")
    if _dab(cfg):
        t1, = st.tabs(["DAB PI control"])
        with t1:
            _dab_pi(cfg, 'th3|pi')
        return
    t1, = st.tabs(["Hysteresis current control"])
    with t1:
        _hysteresis(cfg, 'th3|hys')


#######################################################################################################################
# Theory per mode: mode -> function(cfg)
#######################################################################################################################
THEORY = {0: mode0, 1: mode1, 2: mode2, 3: mode3}
