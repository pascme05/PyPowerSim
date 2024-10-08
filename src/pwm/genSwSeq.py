#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         genSwSeq
# Date:         01.05.2024
# Author:       Dr. Pascal A. Schirmer
# Version:      V.1.0
# Copyright:    Pascal Schirmer
#######################################################################################################################
#######################################################################################################################

#######################################################################################################################
# Function Description
#######################################################################################################################
"""
This function calculates the switching sequence and zero vector for space vector PWM:
Inputs:     1) setup:       includes all simulation variables
Outputs:    1) seq:         switching sequence
            2) k:           zero vectors
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

#######################################################################################################################
# Function
#######################################################################################################################
def genSwSeq(setup):
    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Default SVPWM
    # ==============================================================================
    seq = [[0, 1, 2, 7, 7, 2, 1, 0], [0, 3, 2, 7, 7, 2, 3, 0], [0, 3, 4, 7, 7, 4, 3, 0], [0, 5, 4, 7, 7, 4, 5, 0],
           [0, 5, 6, 7, 7, 6, 5, 0], [0, 1, 6, 7, 7, 6, 1, 0]]
    k = [[0.5, 0.5, 0.5, 0.5, 0.5, 0.5], [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]]

    # ==============================================================================
    # 0127
    # ==============================================================================
    if setup['Par']['PWM']['seq'] == "0127":
        # ------------------------------------------
        # Define Sequence
        # ------------------------------------------
        seq = [[0, 1, 2, 7, 7, 2, 1, 0], [0, 3, 2, 7, 7, 2, 3, 0], [0, 3, 4, 7, 7, 4, 3, 0], [0, 5, 4, 7, 7, 4, 5, 0],
               [0, 5, 6, 7, 7, 6, 5, 0], [0, 1, 6, 7, 7, 6, 1, 0]]

        # ------------------------------------------
        # SVPWM
        # ------------------------------------------
        if setup['Par']['PWM']['zero'] == "SVPWM":
            k = [[0.5, 0.5, 0.5, 0.5, 0.5, 0.5], [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]]

        # ------------------------------------------
        # DPWM0
        # ------------------------------------------
        elif setup['Par']['PWM']['zero'] == "DPWM0":
            k = [[1, 0, 1, 0, 1, 0], [1, 0, 1, 0, 1, 0]]

        # ------------------------------------------
        # DPWM1
        # ------------------------------------------
        elif setup['Par']['PWM']['zero'] == "DPWM1":
            k = [[0, 1, 0, 1, 0, 1], [1, 0, 1, 0, 1, 0]]

        # ------------------------------------------
        # DPWM2
        # ------------------------------------------
        elif setup['Par']['PWM']['zero'] == "DPWM2":
            k = [[0, 1, 0, 1, 0, 1], [0, 1, 0, 1, 0, 1]]

        # ------------------------------------------
        # DPWM3
        # ------------------------------------------
        elif setup['Par']['PWM']['zero'] == "DPWM3":
            k = [[1, 0, 1, 0, 1, 0], [0, 1, 0, 1, 0, 1]]

        # ------------------------------------------
        # DPWMMIN
        # ------------------------------------------
        elif setup['Par']['PWM']['zero'] == "DPWMMIN":
            k = [[1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1]]

        # ------------------------------------------
        # DPWMMAX
        # ------------------------------------------
        elif setup['Par']['PWM']['zero'] == "DPWMMAX":
            k = [[0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]

        # ------------------------------------------
        # Default
        # ------------------------------------------
        else:
            print("WARN: Method not available using SVPWM")
            k = [[0.5, 0.5, 0.5, 0.5, 0.5, 0.5], [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]]

    # ==============================================================================
    # 0121
    # ==============================================================================
    if setup['Par']['PWM']['seq'] == "0121":
        seq = [[0, 1, 2, 1, 1, 2, 1, 0], [0, 3, 2, 3, 3, 2, 3, 0], [0, 3, 4, 3, 3, 4, 3, 0], [0, 5, 4, 5, 5, 4, 5, 0],
               [0, 5, 6, 5, 5, 6, 5, 0], [0, 1, 6, 1, 1, 6, 1, 0]]
        k = [[1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1]]

    # ==============================================================================
    # 7212
    # ==============================================================================
    if setup['Par']['PWM']['seq'] == "7212":
        seq = [[7, 2, 1, 2, 2, 1, 2, 7], [7, 2, 3, 2, 2, 3, 2, 7], [7, 4, 3, 4, 4, 3, 4, 7], [7, 4, 5, 4, 4, 5, 4, 7],
               [7, 6, 5, 6, 6, 5, 6, 7], [7, 6, 1, 6, 6, 1, 6, 7]]
        k = [[0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]

    # ==============================================================================
    # 1012
    # ==============================================================================
    if setup['Par']['PWM']['seq'] == "1012":
        seq = [[1, 0, 1, 2, 2, 1, 0, 1], [3, 0, 3, 2, 2, 3, 0, 3], [3, 0, 3, 4, 4, 3, 0, 3], [5, 0, 5, 4, 4, 5, 0, 5],
               [5, 0, 5, 6, 6, 5, 0, 5], [1, 0, 1, 6, 6, 1, 0, 1]]
        k = [[1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1]]

    # ==============================================================================
    # 2721
    # ==============================================================================
    if setup['Par']['PWM']['seq'] == "2721":
        seq = [[2, 7, 2, 1, 1, 2, 7, 2], [2, 7, 2, 3, 3, 2, 7, 2], [4, 7, 4, 3, 3, 4, 7, 4], [4, 7, 4, 5, 5, 4, 7, 4],
               [6, 7, 6, 5, 5, 6, 7, 6], [6, 7, 6, 1, 1, 6, 7, 6]]
        k = [[0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return [seq, k]
