**Quick Start**
This guide gets you from clone to plots in a few minutes.

**1) Install Dependencies**
1. Create and activate a virtual environment.
2. Install packages from `requirements.txt`.

Example:
```powershell
python -m venv venv
.\venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

**2) Choose a Start Script**
1. Inverter topologies: `start.py`
2. DAB topology: `start_dab.py`

**3) Select a Configuration**
Set `setup['Exp']['conf']` to an Excel file under `config/` without the `.xlsx` extension.
Example in `start.py`:
```python
setup['Exp']['conf'] = "default"
```

**4) Pick the Operating Mode**
Set `setup['Exp']['type']`:
1. `0`: sweep
2. `1`: steady‑state
3. `2`: transient
4. `3`: closed‑loop

Example:
```python
setup['Exp']['type'] = 0
```

**5) Choose the Controlled Output**
Set `setup['Exp']['output']`:
1. `Mi`: modulation index
2. `V`: voltage
3. `I`: current
4. `P`: active power
5. `Q`: reactive power
6. `Phi`: DAB phase shift

**6) Run**
```powershell
python start.py
```
or
```powershell
python start_dab.py
```

**7) Plotting**
Plot control is set by:
1. `setup['Exp']['plot']`
2. `setup['Exp']['plotGen']`
3. `setup['Exp']['save']`

Details are in `docu/PLOTTING.md`.
