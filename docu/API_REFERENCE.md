**API Reference**
This reference lists key modules and functions. It is intended as a map of the codebase, not a full docstring dump.

**Entry Point**
1. `src/main.py`  
   Function `main(setup, path)` executes the full pipeline from load to plot.

**Topology Initialization**
1. `src/topo/initTopo.py`  
   Function `initTopo(name, setup, para)` instantiates the selected topology class.

**Topology Classes**
1. `src/topo/classB2.py`  
2. `src/topo/classB4.py`  
3. `src/topo/classB6.py`  
4. `src/topo/classDAB.py`  

Common methods:
1. `initData()` returns empty time/loss/thermal structures for devices.
2. `initOut()` returns empty time/freq/sweep structures for a run.
3. `out(...)` packages the final output dictionaries.
4. `calcPWM(...)` builds switching sequences.
5. `calcTime(...)` computes time‑domain electrical results.
6. `calcDist(...)` provides analytic distortion (if implemented).
7. `calcCON(...)` updates switching states for closed‑loop control.
8. `initCON()` and `appCON(...)` support controller state and logging.

**Operating Modes**
Inverter:
1. `src/mode/inv/calcSweep.py`
2. `src/mode/inv/calcSteady.py`
3. `src/mode/inv/calcTrans.py`
4. `src/mode/inv/calcClose.py`

DCDC:
1. `src/mode/dcdc/calcSweep_DCDC.py`
2. `src/mode/dcdc/calcSteady_DCDC.py`
3. `src/mode/dcdc/calcTrans_DCDC.py`
4. `src/mode/dcdc/calcClose_DCDC.py`

Each function returns `time`, `freq`, and optionally `sweep`.

**General Utilities**
1. `src/general/genTF.py`  
   Builds transfer functions and state‑space models.
2. `src/general/genLoadInput.py`  
   Computes derived operating points such as `Mi` and `PhiDAB`.
3. `src/general/calcSpec.py` and `src/general/calcFreq.py`  
   Spectral calculations for time‑domain signals.
4. `src/general/calcDistNum.py`  
   Generic numeric distortion metrics for any signal.
5. `src/general/helpFnc.py`  
   Utility functions such as RMS, THD, and `calcDistSignals`.
6. `src/general/sanityCheck.py`  
   Safety checks and bounds.

**Electrical Models**
1. `src/elec/calcElecSwi.py`  
   Switch electrical waveforms.
2. `src/elec/calcLossSwi.py`  
   Switch conduction and switching losses.
3. `src/elec/calcElecCap.py`  
   Capacitor voltage based on ripple current.
4. `src/elec/calcLossCap.py`  
   Capacitor losses from RMS current.
5. `src/elec/calcElecTra.py` and `src/elec/calcLossTra.py`  
   Transformer electrical and loss models.

**Thermal Models**
1. `src/therm/initRC.py`  
   Reduced‑order thermal network parameters.
2. `src/therm/calcTherRC.py`  
   Transient thermal response.

**Plotting**
1. `src/plot/plot.py`  
   Router for generic and topology‑specific plots.
2. `src/plot/gen/*`  
   Generic plots for sweep, steady, transient, closed‑loop.
3. `src/plot/spe/inv/*`  
   Topology‑specific inverter plots.
4. `src/plot/spe/dcdc/*`  
   Topology‑specific DAB plots.

**Result Tables**
1. `src/plot/plotResults.py`  
   Inverter result tables.
2. `src/plot/plotResults_DCDC.py`  
   DCDC result tables.
