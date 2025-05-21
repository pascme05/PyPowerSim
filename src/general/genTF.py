#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         genTF
# Date:         27.04.2024
# Author:       Dr. Pascal A. Schirmer
# Version:      V.1.0
# Copyright:    Pascal Schirmer
#######################################################################################################################
#######################################################################################################################

#######################################################################################################################
# Function Description
#######################################################################################################################
"""
This function checks the parameters for correctness. It also considers general parameters of the machine to assure
smooth calculation.
Inputs:     1) para:    all parameters used in the simulation
            2) setup:   includes all simulation variables
Outputs:    1) out:     the output file includes the transfer functions of the architecture
"""
import numpy as np
#######################################################################################################################
# Import libs
#######################################################################################################################
# ==============================================================================
# Internal
# ==============================================================================

# ==============================================================================
# External
# ==============================================================================
from scipy import signal
import control as ct
import slycot

#######################################################################################################################
# Function
#######################################################################################################################
def genTF(para, setup):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("------------------------------------------")
    print("START: Generate transfer functions")
    print("------------------------------------------")

    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Parameters
    # ==============================================================================
    # ------------------------------------------
    # Load
    # ------------------------------------------
    R = setup['Top']['R']
    L = setup['Top']['L']

    # ------------------------------------------
    # DC-Link
    # ------------------------------------------
    C = para['Cap']['Elec']['con']['C']
    ESR = para['Cap']['Elec']['con']['ESR']

    # ------------------------------------------
    # Transformer
    # ------------------------------------------
    n1 = para['Tra']['Elec']['con']['n1']
    setup['Top']['n1'] = n1
    n2 = para['Tra']['Elec']['con']['n2']
    r1 = para['Tra']['Elec']['con']['r1']
    setup['Top']['r1'] = r1
    r2 = para['Tra']['Elec']['con']['r2']
    setup['Top']['r2'] = r2
    Rc1 = para['Tra']['Elec']['con']['Rc1']
    setup['Top']['Rc1'] = Rc1
    Rc2 = para['Tra']['Elec']['con']['Rc2']
    setup['Top']['Rc2'] = Rc2
    Ll1 = para['Tra']['Elec']['con']['Ll1']
    setup['Top']['Ll1'] = Ll1
    Ll2 = para['Tra']['Elec']['con']['Ll2']
    setup['Top']['Ll2'] = Ll2
    LM = para['Tra']['Elec']['con']['LM']
    M = (n2 / n1) * LM  # Iyer p. 423
    Lm1 = (n1 / n2) * M
    Lm2 = (n2 / n1) * M
    setup['Top']['Ae'] = para['Tra']['Elec']['con']['Ae']

    # ------------------------------------------
    # Output Filter
    # ------------------------------------------
    Rout = setup['Top']['Rout']
    Lout = setup['Top']['Lout']
    Cout = setup['Top']['Cout']

    # ------------------------------------------
    # Input Filter
    # ------------------------------------------
    Rinp = setup['Top']['Rinp']
    Linp = setup['Top']['Linp']
    Cinp = setup['Top']['Cinp']

    # ==============================================================================
    # Variables
    # ==============================================================================
    out = {'TF': {}, 'SS': {}}

    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    # ==============================================================================
    # Input Filter
    # ==============================================================================
    print("INFO: Input filter")
    out['TF']['Inp'] = signal.TransferFunction([ESR*C, 1], [Linp*Cinp, Rinp*Cinp + ESR*C, 1])

    # ==============================================================================
    # DC-Link
    # ==============================================================================
    print("INFO: DC-Link")
    out['TF']['DC'] = signal.TransferFunction([ESR*C, 1], [C, 0])

    # ==============================================================================
    # Output Filter
    # ==============================================================================
    print("INFO: Output filter")
    out['TF']['Out'] = signal.TransferFunction([L, R], [Lout*L*Cout, (Rout*L*Cout + Lout*R*Cout), (Lout + L + R*Rout*Cout), (R + Rout)])

    # ==============================================================================
    # Load
    # ==============================================================================

    # Distinction between the settings for the transformer:
    # Transfer functions for the transformer are based on the simplified equivalent circuit from Iyer, p. 119
    # out['TF']['Load'] is the transfer function that gives i_w1 (primary winding current) from the input voltage v_1
    # other transfer functions are to calculate i_w2 and v_2 for plotting
    if setup['Top']['LD_tra'] == 'NT':
        # No Transformer
        print("INFO: Load")
        out['TF']['Load'] = signal.TransferFunction([1], [L, R])
    elif setup['Top']['LD_tra'] == 'OC':
        # Transformer with open circuit at secondary, no RL load!
        print("INFO: Load is a transformer with open-circuit at secondary")
        out['TF']['Load'] = signal.TransferFunction(np.polymul([-M, 0],[(Ll2 + Lm2), r2+Rc2]), np.polymul([(Ll1 * Lm2 + Ll1 * Ll2 + Lm1 * Ll2 + Lm1 * Lm2 - M * M), (r1 * Ll2 + r1 * Lm2 + Ll1 * r2 + Lm1 * r2 + Ll1*Rc2 + Lm1*Rc2), r1 * r2 + r1*Rc2],[-M, 0]))
        out['TF']['i_w2'] = signal.TransferFunction([-M, 0], [(Ll1 * Lm2 + Ll1 * Ll2 + Lm1 * Ll2 + Lm1 * Lm2 - M * M), (r1 * Ll2 + r1 * Lm2 + Ll1 * r2 + Lm1 * r2 + Ll1*Rc2 + Lm1*Rc2), r1 * r2 + r1*Rc2])
    elif setup['Top']['LD_tra'] == 'SC':
        # Transformer with short circuit at secondary, no RL load!
        print("INFO: Load is a transformer with short-circuit at secondary")
        out['TF']['Load'] = signal.TransferFunction(np.polymul([-M, 0],[(Ll2 + Lm2), +r2]), np.polymul([(Ll1 * Lm2 + Ll1 * Ll2 + Lm1 * Ll2 + Lm1 * Lm2 - M * M), (r1 * Ll2 + r1 * Lm2 + Ll1 * r2 + Lm1 * r2), r1 * r2],[-M, 0]))
        out['TF']['i_w2'] = signal.TransferFunction([-M, 0], [(Ll1 * Lm2 + Ll1 * Ll2 + Lm1 * Ll2 + Lm1 * Lm2 - M * M), (r1 * Ll2 + r1 * Lm2 + Ll1 * r2 + Lm1 * r2), r1 * r2])
    elif setup['Top']['LD_tra'] == 'RL':
        # Transformer with RL+e load connected at the secondary
        print("INFO: Load is a transformer with RL+e load at secondary")
        # TF iw2 = Gi2(v1, e), iw2 is winding current in secondary
        Gi2 = ct.tf([[[L * M, M * (R + Rc2), 0],
                      [-L * M * Rc2 * (Ll1 + Lm1), - L * M * Rc2 * r1 - M * Rc2 * (R + Rc2) * (Ll1 + Lm1),
                       -M * Rc2 * r1 * (R + Rc2), 0]]], [[[L * M * M - (L * Ll2 + L * Lm2) * (Ll1 + Lm1),
                                                            (R + Rc2) * M * M - r1 * (L * Ll2 + L * Lm2) - (
                                                                        Ll1 + Lm1) * (
                                                                        L * r2 + Ll2 * (R + Rc2) + Lm2 * (
                                                                            R + Rc2) + L * Rc2),
                                                            - (R * Rc2 + r2 * (R + Rc2)) * (Ll1 + Lm1) - r1 * (
                                                                        L * r2 + Ll2 * (R + Rc2) + Lm2 * (
                                                                            R + Rc2) + L * Rc2),
                                                            -r1 * (R * Rc2 + r2 * (R + Rc2))],
                                                           [L * M * (L * M * M - (L * Ll2 + L * Lm2) * (Ll1 + Lm1)),
                                                            M * (R + Rc2) * (L * M * M - (L * Ll2 + L * Lm2) * (
                                                                        Ll1 + Lm1)) - L * M * (
                                                                        (- R - Rc2) * M * M + r1 * (
                                                                            L * Ll2 + L * Lm2) + (Ll1 + Lm1) * (
                                                                                    L * r2 + Ll2 * (R + Rc2) + Lm2 * (
                                                                                        R + Rc2) + L * Rc2)),
                                                            - M * (R + Rc2) * ((- R - Rc2) * M * M + r1 * (
                                                                        L * Ll2 + L * Lm2) + (Ll1 + Lm1) * (
                                                                                            L * r2 + Ll2 * (
                                                                                                R + Rc2) + Lm2 * (
                                                                                                        R + Rc2) + L * Rc2)) - L * M * (
                                                                        (R * Rc2 + r2 * (R + Rc2)) * (
                                                                            Ll1 + Lm1) + r1 * (
                                                                                    L * r2 + Ll2 * (R + Rc2) + Lm2 * (
                                                                                        R + Rc2) + L * Rc2)),
                                                            - M * (R + Rc2) * ((R * Rc2 + r2 * (R + Rc2)) * (
                                                                        Ll1 + Lm1) + r1 * (L * r2 + Ll2 * (
                                                                        R + Rc2) + Lm2 * (
                                                                                                       R + Rc2) + L * Rc2)) - L * M * r1 * (
                                                                        R * Rc2 + r2 * (R + Rc2)),
                                                            -M * r1 * (R + Rc2) * (R * Rc2 + r2 * (R + Rc2)), 0]]],
                    inputs=['v1', 'e'], outputs=['iw2'])
        # TF iw1 = Gi2(iw2, e), iw1 is winding current in primary
        Gi1 = ct.tf([[[- L * Ll2 - L * Lm2, - L * r2 - Ll2 * (R + Rc2) - Lm2 * (R + Rc2) - L * Rc2,
                       - R * Rc2 - r2 * (R + Rc2)], [Rc2]]],
                    [[[L * M, M * (R + Rc2), 0], [L * M, M * (R + Rc2), 0]]], inputs=['iw2', 'e'], outputs=['iw1'])
        # Total TF iw1 = G(v1, e)
        out['TF']['Load'] = ct.interconnect([Gi2, Gi1], inplist=['v1', 'e'], outlist=['iw1'])
        out['TF']['i_w2'] = Gi2
        out['TF']['v_2'] = ct.tf([[[-Rc2*L, -Rc2*R], [Rc2]]], [[[L, Rc2+R], [L, Rc2+R]]], inputs=['i_w2', 'e'], outputs=['v_2'])
        out['TF']['i_2'] = ct.tf([[[1], [-1]]], [[[-L, -R], [-L, -R]]], inputs=['v_2', 'e'], outputs=['i_2'])

    ###################################################################################################################
    # Post-Processing
    ###################################################################################################################
    # ==============================================================================
    # Convert-to-State-Space
    # ==============================================================================
    out['SS']['Inp'] = out['TF']['Inp'].to_ss()
    out['SS']['DC'] = out['TF']['DC'].to_ss()
    out['SS']['Out'] = out['TF']['Out'].to_ss()
    if setup['Top']['LD_tra'] != 'RL':
        out['SS']['Load'] = out['TF']['Load'].to_ss()

    
    ###################################################################################################################
    # MSG Out
    ###################################################################################################################
    print("------------------------------------------")
    print("END: Generate transfer functions")
    print("------------------------------------------")

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return out
