**Plotting Guide**
Plotting is controlled from the start scripts using the `setup['Exp']` fields and routed through `src/plot/plot.py`.

**Plot Modes**
1. `plot = 0`: no plots.
2. `plot = 1`: generic plots only.
3. `plot = 2`: generic plots with analytic overlays when available.
4. `plot = 3`: topology‑specific plots in `src/plot/spe/**`.

**Generic Plots**
Generic plots are in `src/plot/gen/`:
1. `plotSweep.py`: modulation or phase‑shift sweep.
2. `plotStat.py`: steady‑state waveforms.
3. `plotTrans.py`: transient waveforms.
4. `plotClose.py`: closed‑loop behavior.

**Topology‑Specific Plots**
1. Inverter plots: `src/plot/spe/inv/`.
2. DAB plots: `src/plot/spe/dcdc/`.
3. DAB closed‑loop plot: `src/plot/spe/dcdc/plotClose_DAB.py`.

**Plotting Extras**
Set `plotGen = 1` to enable:
1. Transfer function Bode plots.
2. Loss model plots.
3. Thermal model plots.

**Result Tables**
Tabular summaries are printed by:
1. `src/plot/plotResults.py` for inverter topologies.
2. `src/plot/plotResults_DCDC.py` for DAB.

**Saving Plots and Results**
Set `save = 1` to store outputs under `results/`. The saving logic lives in `src/general/saveResults.py`.

**Common Plot Controls**
1. `fsim`: affects time resolution and spectrum.
2. `stat.cyc`: controls how many cycles are visible.
3. `plot`: switches between generic and topology‑specific plots.
