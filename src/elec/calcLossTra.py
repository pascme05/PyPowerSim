#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         calcLossTra
# Date:         13.05.2025
# Author:       Max V. Mueller
# Version:      V.1.0
# Copyright:    Pascal Schirmer
#######################################################################################################################
#######################################################################################################################

#######################################################################################################################
# Function Description
#######################################################################################################################
"""
This function calculates the losses of the transformer.
Inputs:     ToDo
Outputs:    ToDo
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
import math
import cmath
import pandas as pd
import numpy as np
from scipy import integrate
from scipy import interpolate
from scipy.special import iv
import matplotlib.pyplot as plt

# Constants:
vacuum_permeability = 1.25663706212e-6  # Vacuum permeability in henries per meter
copper_permeability = 0.999994  # Relative permeability of copper
copper_resistivity = 0.00000001678  # Resistivity of copper in ohm-meters
copper_temperature_coefficient = 0.004041  # Temperature coefficient for copper

#######################################################################################################################
# Helper Functions for round wire
#######################################################################################################################
def modified_bessel_first_kind(order, z):
    return iv(order, z)  # Return the result


def calculate_skin_depth(frequency, temperature):
    """
    Calculates the skin depth for a given frequency and temperature.
    :param frequency: Frequency of the AC current.
    :param temperature: Temperature in degrees Celsius.
    :return: The calculated skin depth.
    """
    # Adjust resistivity based on temperature
    resistivity = copper_resistivity * (1 + copper_temperature_coefficient * (temperature - 20))

    # Calculate skin depth using the formula
    skin_depth = math.sqrt(resistivity / (math.pi * frequency * vacuum_permeability * copper_permeability))
    return skin_depth  # Return the calculated skin depth


def calculate_skin_factor_round_wire(conducting_diameter, outer_diameter, number_conductors, frequency, temperature):
    """
    Calculates the skin factor for a round wire.
    :param conducting_diameter: Diameter of the conducting wire.
    :param outer_diameter: Outer diameter of the wire.
    :param number_conductors: Number of conductors in the wire.
    :param frequency: Frequency of the AC current.
    :param temperature: Temperature in degrees Celsius.
    :return: The calculated skin factor.
    """
    skin_depth = calculate_skin_depth(frequency, temperature)  # Calculate skin depth
    wire_radius = conducting_diameter / 2  # Calculate the radius of the conducting wire
    wire_outer_radius = outer_diameter / 2  # Calculate the outer radius of the wire

    alpha = complex(1, 1)  # Initialize a complex number for calculations
    alpha *= wire_radius / skin_depth  # Scale alpha by the ratio of wire radius to skin depth

    # Calculate the skin factor using the modified Bessel functions
    factor = 0.5 * (alpha * (modified_bessel_first_kind(0, alpha) / modified_bessel_first_kind(1, alpha) + number_conductors * (number_conductors - 1) * (wire_radius ** 2) / (wire_outer_radius ** 2) * modified_bessel_first_kind(1, alpha) / modified_bessel_first_kind(0, alpha))).real

    return factor  # Return the calculated skin factor


#######################################################################################################################
# Helper Function for elliptic integral for rectangular wire
#######################################################################################################################
def comp_ellint_1(x):
    """
    Computes the complete elliptic integral of the first kind. Reference: Zacharias, InductiveDevicesInPE
    x: The parameter for the elliptic integral.
    return: The value of the elliptic integral.
    """

    if x == 0.0:
        return math.pi / 2  # If x is zero, return π/2

    k = abs(x)  # Take the absolute value of x
    m = k * k  # Square of k

    if m == 1.0:
        return float('nan')  # If m is 1, return NaN (not a number)

    a = 1.0  # Initialize a
    g = math.sqrt(1.0 - m)  # Calculate g based on m
    two_n = 1.0  # Initialize two_n

    for i in range(100):  # Iterate up to 100 times
        g_old = g  # Store the old value of g
        a_old = a  # Store the old value of a
        a = 0.5 * (g_old + a_old)  # Update a to the average of g_old and a_old
        g = g_old * a_old  # Update g to the product of g_old and a_old
        two_n += two_n  # Double the value of two_n

        # Break if the change is small enough
        if abs(a_old - g_old) <= (a_old * math.ulp(1.0)):
            break

        g = math.sqrt(g)  # Update g to its square root

    return math.pi / 2 / a  # Return the result of the elliptic integral


#######################################################################################################################
# Function
#######################################################################################################################
def calcLossTra(t, timeAc, T_tr, para, setup):

    ###################################################################################################################
    # Initialisation
    ###################################################################################################################

    # ==============================================================================
    # Variables
    # ==============================================================================
    dt = t[1] - t[0]
    f = np.linspace(0, int(1 / dt), int(len(timeAc['i_1']) / 2))
    fel = setup['Top']['fel']

    # ==============================================================================
    # Output
    # ==============================================================================
    out = pd.DataFrame(columns=['p_L'])

    ###################################################################################################################
    # Pre-Processing
    ###################################################################################################################

    # ==============================================================================
    # Get Fundamental Frequency
    # ==============================================================================
    idx = np.argmin(f - 2 * fel)

    # ==============================================================================
    # Get transformer waveforms
    # ==============================================================================
    i_1 = timeAc['i_1']
    i_w1 = timeAc['i_w1']
    i_2 = timeAc['i_2']
    i_w2 = timeAc['i_w2']
    v_1 = timeAc['v_1']
    v_2 = timeAc['v_2']
    v_w1 = timeAc['v_w1']
    v_w2 = timeAc['v_w2']
    B = timeAc['B']
    dB = timeAc['dB']

    # ==============================================================================
    # RMS winding currents
    # ==============================================================================
    I_w1_rms_sq = np.sum(i_w1 ** 2) / len(i_w1)
    I_w2_rms_sq = np.sum(i_w2 ** 2) / len(i_w2)

    # ==============================================================================
    # Extract Parameters
    # ==============================================================================
    # ------------------------------------------
    # Constant
    # ------------------------------------------
    n1 = para['Tra']['Elec']['con']['n1']
    n2 = para['Tra']['Elec']['con']['n2']
    r1 = para['Tra']['Elec']['con']['r1']
    r2 = para['Tra']['Elec']['con']['r2']
    Rc1 = para['Tra']['Elec']['con']['Rc1']
    setup['Top']['Rc1'] = Rc1
    Rc2 = para['Tra']['Elec']['con']['Rc2']
    setup['Top']['Rc2'] = Rc2
    Ll1 = para['Tra']['Elec']['con']['Ll1']
    Ll2 = para['Tra']['Elec']['con']['Ll2']
    LM = para['Tra']['Elec']['con']['LM']
    M = (n2 / n1) * LM  # Iyer p. 423
    Lm1 = (n1 / n2) * M
    Lm2 = (n2 / n1) * M
    Ae = para['Tra']['Elec']['con']['Ae']
    Ve = para['Tra']['Elec']['con']['Ve']

    # ------------------------------------------
    # Tabular
    # ------------------------------------------

    if setup['Top']['TraMdl'] == "tab":
        # Steinmetz Parameter alpha
        # Matrix
        alpha_2d = interpolate.interp2d(para['Tra']['Elec']['vec']['T_c'].to_numpy(),
                                      para['Tra']['Elec']['vec']['f'].to_numpy(),
                                      para['Tra']['Elec']['tab']['alpha'].to_numpy(), kind='linear')
        # Fundamental Value
        alpha = alpha_2d(T_tr, fel)
        # Static
        alpha = alpha[0]

        # Steinmetz Parameter beta
        # Matrix
        beta_2d = interpolate.interp2d(para['Tra']['Elec']['vec']['T_c'].to_numpy(),
                                        para['Tra']['Elec']['vec']['f'].to_numpy(),
                                        para['Tra']['Elec']['tab']['beta'].to_numpy(), kind='linear')
        # Fundamental Value
        beta = beta_2d(T_tr, fel)
        # Static
        beta = beta[0]

        # Pvsin for calculation of Steinmetz k
        # Matrix
        Pvsin_2d = interpolate.interp2d(para['Tra']['Elec']['vec']['T_c'].to_numpy(),
                                       para['Tra']['Elec']['vec']['f'].to_numpy(),
                                       para['Tra']['Elec']['tab']['Pvsin'].to_numpy(), kind='linear')
        # Fundamental Value
        Pvsin = Pvsin_2d(T_tr, fel)
        # Static
        Pvsin = Pvsin[0] * 1000     # kW --> W
    
    # ------------------------------------------
    # Default
    # ------------------------------------------
    else:
        alpha = para['Tra']['Elec']['con']['alpha']
        beta = para['Tra']['Elec']['con']['beta']
        Pvsin = para['Tra']['Elec']['con']['Pvsin'] * 1000     # kW --> W

    ###################################################################################################################
    # Calculation
    ###################################################################################################################

    # ==============================================================================
    # Winding losses
    # ==============================================================================

    # Calculation for different frequency components, as current is non-sinusoidal:
    n = len(i_w1)
    dT = 1 / setup['Exp']['fsim']

    frequencies = np.fft.fftfreq(n, dT)  # Frequency-vector
    fft_iw1 = np.fft.fft(i_w1)  # FFT
    fft_iw2 = np.fft.fft(i_w2)  # FFT

    half_n = n // 2
    frequencies = frequencies[:half_n]
    magnitude_iw1 = np.abs(fft_iw1[:half_n]) * (2 / n)
    magnitude_iw2 = np.abs(fft_iw2[:half_n]) * (2 / n)

    p_wL1_plusSkin_fft = 0
    p_wL2_plusSkin_fft = 0

    # Differentiate between round and rectangular conductors
    if para['Tra']['Elec']['con']['type_c'] == 'RND':
        # Round wire
        # Skin effect (Zacharias, InductiveDevicesInPE, chapter 9, Albach's model)
        # Parameters
        conducting_diameter = para['Tra']['Elec']['con']['d_c']
        outer_diameter = para['Tra']['Elec']['con']['d_c_o']
        N_cond = para['Tra']['Elec']['con']['N_c']
        for i in range(len(frequencies)):
            Ii_w1_rms_sq = (magnitude_iw1[i] / np.sqrt(2))**2
            Ii_w2_rms_sq = (magnitude_iw2[i] / np.sqrt(2))**2
            fi = frequencies[i]
            if fi == 0:
                Fi_skin = 1
            else:
                Fi_skin = calculate_skin_factor_round_wire(conducting_diameter, outer_diameter, N_cond, fi, T_tr)  # Calculate skin factor
            p_wL1_plusSkin_fft = p_wL1_plusSkin_fft + (r1 * Ii_w1_rms_sq) * Fi_skin
            p_wL2_plusSkin_fft = p_wL2_plusSkin_fft + (r2 * Ii_w2_rms_sq) * Fi_skin
    else:
        # Rectangular wire
        # Skin effect (Zacharias, InductiveDevicesInPE, chapter 9, Kutkut's method)
        # Variables
        conducting_width = para['Tra']['Elec']['con']['w_c']
        conducting_height = para['Tra']['Elec']['con']['h_c']

        b_prime = max(conducting_width, conducting_height) / 2  # Half of the larger dimension
        a_prime = min(conducting_width, conducting_height) / 2  # Half of the smaller dimension

        # Adjust resistivity based on temperature
        resistivity = copper_resistivity * (1 + copper_temperature_coefficient * (T_tr - 20))

        # Calculate frequency parameters
        fl = 3.22 * resistivity / (8 * vacuum_permeability * b_prime * a_prime)  # Low frequency parameter
        fh = (math.pi ** 2 * resistivity) / (4 * vacuum_permeability * (a_prime ** 2)) * (
                comp_ellint_1(math.sqrt(1 - (a_prime ** 2) / (b_prime ** 2))) ** -2)  # High frequency parameter

        # Constants for resistance factor calculation
        gamma = 11
        beta = 5.5
        alpha = 2

        # Basic:
        p_wL = r1*I_w1_rms_sq + r2*I_w2_rms_sq

        # Calculate AC resistance factor for one rms value
        F_skin = (1 + (fel / fl) ** alpha + (fel / fh) ** beta) ** (1 / gamma)

        # Losses
        p_wL_plusSkin = (r1*I_w1_rms_sq + r2*I_w2_rms_sq) + (r1*I_w1_rms_sq + r2*I_w2_rms_sq)*(F_skin - 1)                 # Zacharias p. 86


        for i in range(len(frequencies)):
            Ii_w1_rms_sq = (magnitude_iw1[i] / np.sqrt(2))**2
            Ii_w2_rms_sq = (magnitude_iw2[i] / np.sqrt(2))**2
            fi = frequencies[i]
            Fi_skin = (1 + (fi / fl) ** alpha + (fi / fh) ** beta) ** (1 / gamma)
            p_wL1_plusSkin_fft = p_wL1_plusSkin_fft + (r1 * Ii_w1_rms_sq) * Fi_skin
            p_wL2_plusSkin_fft = p_wL2_plusSkin_fft + (r2 * Ii_w2_rms_sq) * Fi_skin

    # Proximity losses:
    p_w_prox = 0        # not considered so far, as strongly dependent on geometry, can be implemented via tabular ac-resistance in the future

    # Total winding losses
    p_wL_total_fft = p_wL1_plusSkin_fft + p_wL2_plusSkin_fft + p_w_prox

    # ==============================================================================
    # Core losses
    # ==============================================================================
    # Modified Steinmetz-Equation (MSE)
    B_pk = np.max(B)
    B_pk2pk = np.max(B) - np.min(B)

    f_eq = (2/(B_pk2pk**2*np.pi**2)) * integrate.trapezoid(dB**2, dx=dt)
    k = Pvsin/(f_eq**(alpha-1)*B_pk**beta*fel)
    p_cL = (k*f_eq**(alpha-1)*B_pk**beta*fel)*Ve

    # Create vectors of constant power loss
    p_wL1_plusSkin_fft = p_wL1_plusSkin_fft * np.ones(np.shape(i_w1))
    p_wL2_plusSkin_fft = p_wL2_plusSkin_fft * np.ones(np.shape(i_w1))
    p_cL = p_cL * np.ones(np.shape(i_w1))

    # Output
    out['p_wL_1'] = p_wL1_plusSkin_fft          # Winding losses primary
    out['p_wL_2'] = p_wL2_plusSkin_fft          # Winding losses secondary
    out['p_cL'] = p_cL                          # Core losses
    out['p_L_total'] = p_wL_total_fft + p_cL    # Total losses

    ###################################################################################################################
    # Post-Processing
    ###################################################################################################################

    ###################################################################################################################
    # Return
    ###################################################################################################################
    return out
