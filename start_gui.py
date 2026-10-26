#######################################################################################################################
#######################################################################################################################
# Title:        Power Electronics Simulation Toolkit (PyPowerSim)
# Topic:        Power Electronics
# File:         start_gui.py
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
This function starts the browser GUI of the toolkit locally (http://localhost:8501). All parameters of the operating
modes are edited in the GUI and stored as JSON under config/; the input data are read from the workbooks. Requires the
packages in gui/requirements.txt.
"""

#######################################################################################################################
# Import external libs
#######################################################################################################################
import os
import subprocess
import sys

#######################################################################################################################
# Start
#######################################################################################################################
if __name__ == '__main__':
    root = os.path.dirname(os.path.abspath(__file__))
    print("GUI: http://localhost:8501")
    try:
        import streamlit  # noqa: F401
    except ImportError:
        sys.exit("ERROR: streamlit is not installed, run: pip install -r gui/requirements.txt")
    if not os.path.isfile(os.path.join(root, 'gui', 'schema.json')):
        subprocess.check_call([sys.executable, os.path.join(root, 'gui', 'build_schema.py')], cwd=root)
    sys.exit(subprocess.call([sys.executable, '-m', 'streamlit', 'run', os.path.join(root, 'gui', 'app.py')]
                             + sys.argv[1:], cwd=root))
