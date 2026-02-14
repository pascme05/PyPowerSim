**Tutorials**
These short tutorials walk through typical workflows for inverter and DAB analyses.

**Tutorial 1: Inverter Sweep (B4)**
1. Open `start.py`.
2. Set `setup['Exp']['conf'] = "default"` or another inverter config.
3. Set `setup['Top']['sourceType'] = "B4"` in the config file.
4. Set `setup['Exp']['type'] = 0` for sweep.
5. Run `python start.py`.
6. Inspect sweep plots for modulation‑dependent distortion.

**Tutorial 2: Inverter Steady‑State (B6)**
1. Open `start.py`.
2. Set `setup['Exp']['conf'] = "default"` with a B6 topology.
3. Set `setup['Exp']['type'] = 1`.
4. Run `python start.py`.
5. Use `plotResults` output tables for switch losses and temperatures.

**Tutorial 3: Inverter Transient**
1. Open `start.py`.
2. Set `setup['Exp']['type'] = 2`.
3. In the config, set `setup['Dat']['trans']['tmax']` to a meaningful time window.
4. Run and inspect thermal evolution plots.

**Tutorial 4: DAB Sweep**
1. Open `start_dab.py`.
2. Set `setup['Exp']['conf'] = "default_DCDC"`.
3. Set `setup['Exp']['type'] = 0`.
4. Run `python start_dab.py`.
5. Inspect DAB sweep plots for phase‑shift distortion and power.

**Tutorial 5: DAB Steady‑State**
1. Open `start_dab.py`.
2. Set `setup['Exp']['type'] = 1`.
3. Adjust `PhiDAB` in the config if you want a different operating point.
4. Run and inspect loss tables in `plotResults_DCDC`.

**Tutorial 6: Parameter Variation**
1. Duplicate a config file in `config/`.
2. Adjust one parameter at a time, such as `fs`, `PhiDAB`, or `Vdc`.
3. Re‑run the same analysis type and compare plots.
