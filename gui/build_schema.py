#######################################################################################################################
#######################################################################################################################
# Title:        Power Electronics Simulation Toolkit (PyPowerSim)
# Topic:        Power Electronics
# File:         build_schema.py
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
This function converts the config workbook (adapter.CONFIG_WORKBOOK) into the parameter schema of the GUI
(gui/schema.json) and the default configuration (config/default.json). The schema holds name, description, unit, type
and selectable options of every parameter, the GUI forms are generated from it. Only needed if parameters are added to
the workbook, the GUI itself does not read the workbook.
Usage:      python gui/build_schema.py [path/to/config.xlsx]
"""

#######################################################################################################################
# Import libs
#######################################################################################################################
import json
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from gui import adapter as A


#######################################################################################################################
# Helper
#######################################################################################################################
def _clean(val):
    """Excel cell to JSON value (text with leading zeros such as '0127' stays text)."""
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return None
    if isinstance(val, str):
        val = val.strip()
        if re.fullmatch(r'0\d+', val):
            return val
        try:
            val = float(val)
        except ValueError:
            return val
    if isinstance(val, (bool, np.bool_)):
        return bool(val)
    if isinstance(val, (int, np.integer)):
        return int(val)
    val = float(val)
    return int(val) if val.is_integer() and abs(val) < 1e15 else val


MARK = re.compile(r'(?:^|(?<=\n)|(?<=[,;:|] ))(?:\((-?[\w.]+)\):|(-?\d+|[A-Z][a-z]{2})\))\s')


def _options(desc):
    """
    Selectable options from a description: '1) Normal\n2) Uniform', '(0): off,\n(1): on', '(SV): space vector' or
    '(0: off, 1: on)'. Returns the remaining description and the options {key: label} (None if no options).
    """
    m = re.search(r'\((-?\d+):\s*([^,()]+),\s*(-?\d+):\s*([^()]+)\)', desc)
    if m:
        opts = {m.group(1): m.group(2).strip(), m.group(3): m.group(4).strip()}
        return (desc[:m.start()] + desc[m.end():]).strip(), opts

    marks = list(MARK.finditer(desc))
    if len(marks) < 2:
        return desc, None
    opts = {}
    for i, mk in enumerate(marks):
        end = marks[i + 1].start() if i + 1 < len(marks) else len(desc)
        label = desc[mk.end():end].split('\n\n')[0]                    # a blank line ends the option (heading)
        opts[mk.group(1) or mk.group(2)] = label.strip().rstrip(',;|').strip()
    head = desc[:marks[0].start()].strip().rstrip(':').strip()
    return head, opts


#######################################################################################################################
# Function
#######################################################################################################################
def build_schema(filename):
    ###################################################################################################################
    # Loading
    ###################################################################################################################
    raw = pd.read_excel(filename, sheet_name=None)

    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    schema, default = {}, {}
    for sheet, (section, nested) in A.SHEETS.items():
        params = []
        for _, r in raw[sheet].iterrows():
            var, desc = _clean(r['Variable']), _clean(r['Description'])

            # Continuation rows (description only) are appended to the previous parameter
            if var is None:
                if desc is not None and params:
                    params[-1]['desc'] += '\n' + str(desc)
                continue

            var = str(var)
            cat = str(_clean(r['Category']) or section)
            val, unit = _clean(r['Value']), _clean(r['Unit'])
            unit = None if unit in (None, '-') else str(unit)
            head, opts = _options(str(desc or ''))
            if (section, cat, var) in A.NO_OPTIONS:
                head, opts = str(desc or ''), None
            extra = getattr(A, 'OPTIONS', {}).get((section, cat, var))
            if extra:                                                  # options used by the code, missing in the text
                opts = {**(opts or {}), **{str(k): v for k, v in extra.items()}}
                opts = {k: v for k, v in opts.items() if v is not None} or None    # None removes a wrong option

            # Type (digit keys of a text value such as '0127' stay text)
            text = isinstance(val, str) or (section, cat, var) in A.STRS
            if opts and not text and all(re.fullmatch(r'-?\d+', k) for k in opts):
                typ, opts = 'int', {int(k): v for k, v in opts.items()}
            elif opts:
                typ = 'str'
            elif isinstance(val, bool):
                typ = 'bool'
            elif (section, cat, var) in A.INTS:
                typ = 'int'
            elif text:
                typ = 'str'
            else:
                typ = 'float'
            if typ == 'str' and opts and val is not None and not isinstance(val, str):
                val = str(val)
            name = getattr(A, 'LABELS', {}).get((section, cat, var)) or _clean(r['Name']) or var
            percent = unit == '%' and getattr(A, 'PERCENT_AS_FRACTION', True)

            params.append({'cat': cat, 'var': var, 'name': name, 'desc': head, 'unit': unit,
                           'type': typ, 'percent': percent, 'options': opts, 'nested': nested,
                           'optional': (section, cat, var) in A.OPTIONAL, 'default': val})
            if nested:
                default.setdefault(section, {}).setdefault(cat, {})[var] = val
            else:
                default.setdefault(section, {})[var] = val
        schema.setdefault(section, []).extend(params)

    ###################################################################################################################
    # Post-processing
    ###################################################################################################################
    default[A.RUN_SECTION] = {**A.RUN, **default.get(A.RUN_SECTION, {})}

    return schema, default


#######################################################################################################################
# Main
#######################################################################################################################
if __name__ == '__main__':
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(A.ROOT, A.CONFIG_WORKBOOK)
    schema, default = build_schema(src)
    with open(os.path.join(A.ROOT, 'gui', 'schema.json'), 'w', encoding='utf-8') as f:
        json.dump(schema, f, indent=1, ensure_ascii=False)
    os.makedirs(os.path.join(A.ROOT, A.CONFIG_DIR), exist_ok=True)
    with open(os.path.join(A.ROOT, A.CONFIG_DIR, 'default.json'), 'w', encoding='utf-8') as f:
        json.dump(default, f, indent=2, ensure_ascii=False)

    # Further ready-to-load configs from other workbooks of the same layout (adapter.EXTRA_CONFIGS)
    for name, (wb, run) in getattr(A, 'EXTRA_CONFIGS', {}).items():
        _, cfg = build_schema(os.path.join(A.ROOT, wb))
        cfg[A.RUN_SECTION].update(run)
        with open(os.path.join(A.ROOT, A.CONFIG_DIR, name + '.json'), 'w', encoding='utf-8') as f:
            json.dump(cfg, f, indent=2, ensure_ascii=False)
        print(f"INFO: config {name}.json written from {wb}")
    for sec, ps in schema.items():
        n_opt = sum(1 for p in ps if p['options'])
        print(f"INFO: {sec:<14} {len(ps):>4} parameters ({n_opt} with options)")
    print(f"DONE: schema with {sum(len(v) for v in schema.values())} parameters written")
