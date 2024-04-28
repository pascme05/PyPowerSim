#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         template
# Date:         27.04.2024
# Author:       Dr. Pascal A. Schirmer
# Version:      V.0.2
# Copyright:    Pascal Schirmer
#######################################################################################################################
#######################################################################################################################

#######################################################################################################################
# Function Description
#######################################################################################################################
"""
This function checks the parameters for correctness. It also considers general parameters of the machine to assure
smooth calculation.
Inputs:     1) setup:   includes all simulation variables
Outputs:    1) setup:   updated setup file
"""

#######################################################################################################################
# Import libs
#######################################################################################################################
# ==============================================================================
# Internal
# ==============================================================================

# ==============================================================================
# External
# ==============================================================================
import numpy as np
import math
import cmath


#######################################################################################################################
# Function
#######################################################################################################################
def genLoadInput(setup):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("------------------------------------------")
    print("START: Generate Load Input")
    print("------------------------------------------")

    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Parameters
    # ==============================================================================
    # ------------------------------------------
    # Input/Output
    # ------------------------------------------
    Mi = setup['Dat']['stat']['Mi']
    Vdc = setup['Dat']['stat']['Vdc']
    phiO = math.radians(setup['Dat']['stat']['phi'])

    # ------------------------------------------
    # Load
    # ------------------------------------------
    R = setup['Top']['R']
    L = setup['Top']['L']
    E = setup['Top']['E']
    phiE = math.radians(setup['Top']['phiE'])

    # ------------------------------------------
    # Operating point
    # ------------------------------------------
    fel = setup['Top']['fel']

    # ==============================================================================
    # Topology
    # ==============================================================================
    if setup['Top']['sourceType'] == "B2":
        paraV = 2
    elif setup['Top']['sourceType'] == "B4":
        paraV = 1
    elif setup['Top']['sourceType'] == "B6":
        paraV = 2
    else:
        paraV = 2

    ###################################################################################################################
    # Pre-Processing
    ###################################################################################################################
    # ==============================================================================
    # Load Parameters
    # ==============================================================================
    magZ = np.sqrt(R ** 2 + (2 * np.pi * fel * L) ** 2)
    angZ = math.atan2(2 * np.pi * fel * L, R)
    Z = complex(R, 2 * np.pi * fel * L)
    EMF = complex(E * np.cos(phiE), E * np.sin(phiE)) * Mi

    # ==============================================================================
    # Printing Load
    # ==============================================================================
    print(f"INFO: Complex load Z with magnitude {magZ:.2f} (V/A) and angle {math.degrees(angZ):.2f} (deg)")
    print(f"INFO: Back EMF with magnitude {E:.2f} (V/A) and angle {math.degrees(phiE):.2f} (deg)")

    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    # ==============================================================================
    # Modulation index is controlled
    # ==============================================================================
    if setup['Exp']['output'] == 'Mi':
        print(f"INFO: Modulation index controlled mode with Mi {setup['Dat']['stat']['Mi']:.2f} (p.u.)")

    # ==============================================================================
    # Voltage is controlled
    # ==============================================================================
    elif setup['Exp']['output'] == 'V':
        print(f"INFO: Voltage controlled mode with Vo {setup['Dat']['stat']['Vo']:.2f} (V)")
        Mi = setup['Dat']['stat']['Vo'] / setup['Dat']['stat']['Vdc'] * paraV

    # ==============================================================================
    # Current is controlled
    # ==============================================================================
    elif setup['Exp']['output'] == 'I':
        print(f"INFO: Current controlled mode with Io {setup['Dat']['stat']['Io']:.2f} (A)")
        Vo = abs(setup['Dat']['stat']['Io'] * Z * np.sqrt(2) + EMF)
        Mi = Vo / setup['Dat']['stat']['Vdc'] * paraV

    # ==============================================================================
    # Active power is controlled
    # ==============================================================================
    elif setup['Exp']['output'] == 'P':
        print(f"INFO: Active power controlled mode with Po {setup['Dat']['stat']['Po']:.2f} (W)")
        Io = np.sqrt(setup['Dat']['stat']['Po'] / R) * np.sqrt(2)
        Io = complex(Io * np.cos(angZ), Io * np.sin(angZ))
        Mi = abs(Io * Z + EMF) / setup['Dat']['stat']['Vdc'] * paraV

    # ==============================================================================
    # Reactive power is controlled
    # ==============================================================================
    elif setup['Exp']['output'] == 'Q':
        print(f"INFO: Reactive power controlled mode with Qo {setup['Dat']['stat']['Qo']:.2f} (W)")
        Io = np.sqrt(setup['Dat']['stat']['Qo'] / (2 * np.pi * fel * L)) * np.sqrt(2)
        Io = complex(Io * np.cos(angZ), Io * np.sin(angZ))
        Mi = abs(Io * Z + EMF) / setup['Dat']['stat']['Vdc'] * paraV

    # ==============================================================================
    # Default
    # ==============================================================================
    else:
        print(f"INFO: Modulation index controlled mode with Mi {setup['Dat']['stat']['Mi']:.2f} (p.u.)")

    ###################################################################################################################
    # Post-Processing
    ###################################################################################################################
    # ==============================================================================
    # Output Values
    # ==============================================================================
    Vo = Mi * Vdc / paraV
    Io = (complex(Vo, phiO) - EMF) / Z / np.sqrt(2)
    So = Vo * Io / np.sqrt(2)
    Po = So.real
    Qo = So.imag

    # ==============================================================================
    # Printing
    # ==============================================================================
    print(f"INFO: Modulation index being equal to {Mi:.2f} (p.u.)")
    print(f"INFO: Output rms voltage being equal to {Vo:.2f} (V) and {phiO:.2f} (deg)")
    print(
        f"INFO: Output rms current being equal to {cmath.polar(Io)[0]:.2f} (A) and {math.degrees(cmath.polar(Io)[1]):.2f} (deg)")
    print(f"INFO: Output power being equal to {Po:.2f} (W) and {Qo:.2f} (VAr)")

    # ==============================================================================
    # Output
    # ==============================================================================
    setup['Dat']['stat']['Vo'] = Vo
    setup['Dat']['stat']['Io'] = Io
    setup['Dat']['stat']['Po'] = Po
    setup['Dat']['stat']['Qo'] = Qo
    setup['Dat']['stat']['Mi'] = Mi

    ###################################################################################################################
    # MSG Out
    ###################################################################################################################
    print("------------------------------------------")
    print("END: Generate Load Input")
    print("------------------------------------------")

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return setup
