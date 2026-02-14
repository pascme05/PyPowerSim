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

**What To Look For (Examples)**
1. Modulation / switching plots show `s(t)` waveforms and their spectra.
2. Current plots reveal RMS values and harmonic content.
3. Voltage plots indicate ripple and harmonics on the DC link and AC side.
4. Loss and thermal plots show device stress and steady‑state operating points.

Example figures are already included in `docu/images` and referenced in `README.md`, such as:
1. `docu/images/img.png` (modulation function)
2. `docu/images/img_1.png` (current plots)
3. `docu/images/img_2.png` (voltage plots)
4. `docu/images/steady1.png`, `docu/images/steady2.png` (steady‑state results)

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

**DAB Plot Notes**
1. Primary and secondary results are often shown together for direct comparison.
2. Distortion plots use `PhiDAB` on the x‑axis instead of modulation index.
3. DC‑side plots show both `v_dc_pri` and `v_dc_sec` when available.

**See Also**
1. Inline plot walkthroughs with explanations are in `docu/EXAMPLE_RESULTS.md`.
