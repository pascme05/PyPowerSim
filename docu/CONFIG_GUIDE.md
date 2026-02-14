**Configuration Guide**
PyPowerSim uses Excel configuration files located in `config/`. These files define all experiment, data, topology, and parameter settings. The file is read by `src/data/loadSetup.py`.

**Sheets**
1. `Exp`: experiment‑level settings such as plot and save options.
2. `Dat`: analysis‑level data for steady‑state and transient runs.
3. `Top`: topology parameters such as load and source configuration.
4. `Par`: PWM, control, electrical, and thermal model settings.

**Common Columns**
1. `Name`: descriptive name of the variable.
2. `Category`: group such as `stat` or `trans`.
3. `Description`: explanation and valid options.
4. `Variable`: name used in code.
5. `Value`: editable value used in simulation.
6. `Unit`: unit of the variable.

**Exp Sheet Essentials**
1. `output`: controlled output, such as `Mi`, `V`, `I`, `P`, `Q`, or `Phi`.
2. `type`: operating mode, `0` to `3` for sweep, steady, transient, closed‑loop.
3. `plot`: plot selection, `0` to `3` for none, generic, analytic, or topology‑specific.
4. `plotGen`: generic plots toggle.
5. `save`: results saving toggle.
6. `fsim`: simulation sampling frequency.
7. `eps`: numerical epsilon for timing comparisons.
8. `therFeed`: thermal feedback enable.

**Dat Sheet Essentials**
1. `stat.cyc`: number of fundamental cycles for steady‑state/sweep.
2. `stat.W`: number of sweep points.
3. `stat.Mi`: modulation index.
4. `stat.PhiDAB`: DAB phase shift (deg or rad depending on config source).
5. `stat.Vdc`: DC‑link voltage.
6. `stat.Tc`, `stat.Tj`: temperatures for steady‑state.
7. `trans.tmax`: transient simulation duration.
8. `trans.Tc`, `trans.Tj`: initial temperatures for transient runs.

**Top Sheet Essentials**
1. `sourceType`: `B2`, `B4`, `B6`, or `DAB`.
2. `fel`: electrical output frequency.
3. `R`, `L`: load parameters.
4. `E`, `phiE`: back‑EMF settings for inverter mode.
5. `Rout`, `Lout`, `Cout`: output filter parameters if used.
6. `Rinp`, `Linp`, `Cinp`: input filter parameters if used.

**Par Sheet Essentials**
1. `Par.PWM`: switching settings such as `fs`, `type`, `td`, `tmin`.
2. `Par.Cont`: controller settings, such as `type`, `fc`, `Kp`, `Ki`, `hys`.
3. `Par.Elec`: electrical model settings for switch and capacitor.
4. `Par.Ther`: thermal model settings and coupling.

**DAB Notes**
1. The DAB operating point is driven by `PhiDAB`.
2. `output = Phi` is recommended when you want direct phase‑shift control.
3. `N` and `Lk` are defined in the transformer parameter file referenced by `Exp.Trafo`.

**Tips**
1. Keep `stat.cyc` at 4 or more for clean spectral results.
2. Ensure `fs` is an integer multiple of `fel` for synchronous PWM.
3. Use `sanityCheck` warnings to adjust simulation stability and memory use.
