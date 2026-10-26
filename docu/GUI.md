# PyPowerSim GUI — User Guide

The PyPowerSim GUI is a local web application (Streamlit) for the toolkit. All parameters of the config workbook and
the settings of `start_INV.py` / `start_DCDC.py` are edited in the browser and stored as JSON; the device parameters are
still read from the workbooks in `para/`. The GUI runs on your own machine only (`http://localhost:8501`, not reachable
from the network). The start scripts keep working unchanged.

## 1. Installation and start

```
pip install -r requirements.txt          # toolkit
pip install -r gui/requirements.txt      # GUI: streamlit, plotly, openpyxl
python start_gui.py                      # opens http://localhost:8501
```

`start_gui.py` starts `streamlit run gui/app.py` with the repository root as working directory. Stop the GUI with
`Ctrl+C` in the terminal. Server settings (port, localhost only, theme) are in `.streamlit/config.toml`.

## 2. Layout

| Area | Content |
|---|---|
| **Sidebar** | Expert mode switch, configuration manager (load / save / download / import / reset), runs of this session, navigation |
| **Setup** | Run settings of the start scripts (name, controlled quantity, plots, saving), topology and the device files grouped into *Inverter* and *DC-DC converter*, simulation settings (sheet *Exp*) |
| **Input data** | Device workbooks of the inverter or the DC-DC converter (the group used by the topology first) with charts and editor, plus a browser of all workbooks in `para/Swi`, `para/Cap`, `para/Tra` |
| **0 Sweep · 1 Steady state · 2 Transient · 3 Closed loop** | One page per operating mode (`setup['Exp']['type']`) with the views *Workflow*, *Parameters*, *Inputs*, *Theory*, *Run*, *Results* |
| **Results browser** | All folders in `results/`, incl. runs of the start scripts |
| **Help** | Tutorials, this guide and the parameter reference |

### 2.1 Expert mode

**Off (default):** only the main parameters are shown, and only those that have an effect with the current settings,
e.g. the set point *Vo* only if the controlled quantity is *V*, *PhiDAB* only for the DAB, the filter elements only if
the filter is switched on, the controller only in the closed-loop mode, the transient settings only in modes 2 and 3.
**On:** all parameters, inactive ones are marked *· inactive*. Hidden parameters keep their values and are used by the
simulation. Help texts, basic parameters and rules are in `gui/core/help.py`.

### 2.2 Inverter and DC-DC converter

The topology (*Parameters › Topology & load › Topology*) selects B2, B4, B6 or DAB. Two configurations are prepared:

| Config | Origin | Content |
|---|---|---|
| `default` | `config/default.xlsx` + `start_INV.py` | B6 inverter, SVPWM, sweep mode |
| `default_DCDC` | `config/default_DCDC.xlsx` + `start_DCDC.py` | DAB, single phase shift, steady state, fsim = 100 MHz |

Load them in the sidebar (*Saved configs*). On the Setup page the device files are shown in two boxes: **Inverter**
(switches, DC-link capacitor; used by B2/B4/B6) and **DC-DC converter** (switches of the primary and secondary bridge,
input and output capacitor, transformer; used by the DAB). A badge shows which group the current topology uses; the
other group keeps its files but is not read.

### 2.3 Device workbooks: charts and editor

*Input data* and the *Inputs* view of the modes show every device workbook in four tabs:

| Tab | Content |
|---|---|
| **Parameters** | Ratings of the readme sheet and the parameter table of the sheet *electrical* (symbol, typical value, unit, model, vector Value-1 … 10, description) |
| **Characteristics** | The look-up tables as curves over the x-axis vector, one curve per junction temperature: switch v_CE, E_on, E_off, E_rec, v_F, capacitances; capacitor C, tan δ, ripple factor, ESR; transformer winding resistances, magnetising inductance, core losses, stray capacitances. Empty tables are listed |
| **Thermal network** | Zth(t) of the Foster networks (junction → case of transistor and diode, case → ambient …) with a table of Rth [K/W], Cth [Ws/K] and τ [s] |
| **Edit & save** | Editable parameter table, look-up tables and thermal table. *Save* writes `para/<folder>/<name>_edit.xlsx` (or another name; overwriting needs a confirmation), loads it with the loader of the toolkit as a check and offers *Use in the configuration*. *Discard* resets the edits |

The edits are shown in the other tabs at once (badge *unsaved edits*). The saved workbook keeps the layout of the
original; formulas of the template (e.g. axis vectors `=F$16`) are written as values because the loaders read the
cached cell values.

## 3. Configurations

The configuration is the nested dictionary `setup` of the toolkit (sections *Exp, Dat, Top, Par*). *Save* writes
`config/<name>.json`, *Import* reads a `.json` or any config workbook (`config/*.xlsx`, e.g. the journal or compare
setups). For a run the GUI writes the configuration into `config/_gui_run.xlsx` (a copy of `config/default.xlsx` with
the GUI values), sets `setup['Exp']['conf'] = '_gui_run'` and calls `main()` — the toolkit loads it with its own loader.

After changing the config workbook template regenerate the schema and the default configs:

```
python gui/build_schema.py
```

Options that the code accepts but the workbook does not list (e.g. the *DAB* topology, the B6 modulation methods and
switching sequences) are added in `gui/adapter.py` (`OPTIONS`).

## 4. Theory and tutorials

The view **Theory** of every mode page explains the models with equations and interactive examples. The explorers
start from the current configuration but do not change it; where possible they call the toolkit itself.

| Mode | Explorers |
|---|---|
| 0 Sweep | **PWM explorer**: one period of modulation function, carrier, leg voltage and current computed with `calcRef`/`calcPWM`/`calcTime` for any Mi, PWM method and pulse number, with spectrum and THD; **Fundamental over Mi** (linear range, overmodulation, six-step limit). DAB: power transfer over the phase shift with current waveform and ZVS |
| 1 Steady state | **Losses of a switch**: conduction and switching losses averaged over a sinusoidal current from the device characteristics; **Temperature ripple** of the junction with `calcTherRC` over the fundamental frequency. DAB: power transfer and steady temperature |
| 2 Transient | **Heating of the devices** after a loss step (Zth, time constants, settling at tmax); **Electrical start-up** of the RL load |
| 3 Closed loop | **Hysteresis current control** of one leg with `conHys` (band, fc, switching frequency, error). DAB: PI control of the phase shift |

Note: in B2, B4 and B6 the controller setting *PI* also runs the hysteresis control (`calcCON` of the topology
classes); the PI controller (Kp, Ki, saturation) is implemented for the DAB only.

**Help › Tutorials** contains guided scenarios (configs in `config/tutorials/`): sweep of the B6 inverter, losses and
temperatures of the IGBT inverter, discontinuous PWM, DAB power transfer and losses, heating after switch-on and
hysteresis current control. *Load tutorial* replaces the current configuration and opens the mode page; the tab *What
you should see* lists the results of a run of the tutorial config.

## 5. Running a simulation

The *Run* view shows a summary and the button **Run**. The console output is streamed. Matplotlib figures of the toolkit
are captured (no windows open) and shown under *Results › Figures & files*. With *Save results = 1* the toolkit writes
`results/<Name>/`; the GUI adds `gui_run.pkl` / `gui_run.json` so the run can be reopened in the results browser
(signals are thinned to 40 000 samples for storage).

Run times with the templates: steady state about 1 min, sweep about 40 s, DAB steady state about 45 s, transient
about 2–3 min. The closed loop simulates tmax × fc control periods one by one (template: 0.5 s × 20 kHz = 10 000
periods, more than 40 min); for a first look set *Parameters › Transient › tmax* to e.g. 0.05 s (about 3 min).

## 6. Results

| View | Content |
|---|---|
| **Key figures** | Output power, total losses, efficiency (definitions as in the results tables of the toolkit), maximum temperature, THD of the phase current (inverter, last period) |
| **Waveforms** | Switching functions, AC and DC side signals; evaluation window *Last period* (stationary) or *Full simulation*. In the transient mode the toolkit stores the AC/DC signals of the last period only, losses and temperatures over tmax |
| **Losses** | Mean conduction and switching losses of transistor and diode per switch, capacitor and transformer losses |
| **Temperatures** | Temperatures of transistors, diodes, cases, capacitors and transformer with mean, maximum and swing |
| **Spectrum** | Amplitude spectrum over the harmonic order, logarithmic scale |
| **Operating point** | The setup values after the toolkit derived the operating point (Po, Qo, Io, …) |
| **Sweep** (mode 0) | RMS values, fundamentals and distortion over the swept parameter, numerical vs analytical. The quantities `*_thd` of the toolkit are the RMS value of the harmonics in V or A (`src/general/calcDistNum.py`), not a ratio; the ratio is X_thd / X_v1_eff |

The tables list units for every column; waveforms come with a table of mean, RMS, minimum and maximum per signal and
the spectrum with the largest components (order, frequency, amplitude, share of the fundamental).

Every chart has a line *How to read this chart* (texts in `gui/core/chart_help.py`). The losses in the GUI are the sum
of **all** components; the DAB console table of `src/plot/plotResults_DCDC.py` currently omits S4 and S8 (slices
`P_ALL[0:3]`, `P_ALL[4:7]`), so its total losses are lower and its efficiency higher than in the GUI.

## 7. Files of the GUI

| File | Purpose |
|---|---|
| `start_gui.py` | Starts the GUI |
| `gui/adapter.py` | Connection to the toolkit: sheets, run settings, modes, device folders, call of `main()` (records the outputs of the `calc*` functions) |
| `gui/build_schema.py` | Generates `gui/schema.json`, `config/default.json` and `config/default_DCDC.json` from the config workbooks |
| `gui/core/help.py` | Help texts, basic parameters, activity rules |
| `gui/core/viz_results.py` | Interactive results |
| `gui/core/workflow.py` | Workflow charts |
| `gui/core/chart_help.py` | Explanations of the charts |
| `gui/core/devices.py`, `gui/core/viz_inputs.py` | Reading, charts, editing and saving of the device workbooks |
| `gui/core/theory.py`, `gui/core/toolkit.py` | Theory explorers; pre-processing of the toolkit for one period |
| `gui/core/tutorials.py`, `gui/views/tutorials.py` | Tutorials (`python -m gui.core.tutorials` rewrites `config/tutorials/*.json`) |
| `gui/core/*.py`, `gui/views/*.py` | Generic GUI (forms, runner, files, style, pages) |

## 8. Troubleshooting

| Problem | Remedy |
|---|---|
| Port 8501 in use | Stop the other instance or change `server.port` in `.streamlit/config.toml` |
| Changes in the code are not visible | Restart `start_gui.py` |
| Closed loop stopped with "Length of values (10000) does not match length of index (10001)" | Fixed in `src/mode/inv/calcClose.py`: switching functions shorter than one fundamental period are padded before `out()` |
| A saved device workbook is rejected by the loader | The check after *Save* shows the console output and the error of the loader; correct the named cells in *Edit & save* and save again |
| A run seems to hang | Check fsim × cyc (steady state) or fsim × tmax (transient): memory grows with the number of samples; the console log shows the sanity check warnings |
| A config from an older version misses parameters | Missing parameters are filled with the defaults on loading |
| Interactive charts missing for a run | Runs of the start scripts have no stored data; their figures and files are shown |
