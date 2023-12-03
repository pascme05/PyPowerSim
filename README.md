# PyPowerSim
PyPowerSim is a simple python toolkit for evaluation of standard power converter topologies.
The current version includes simulation architectures for a half-bridge (B2), a 
full-bridge (B4), and a three-phase full-bridge converter (B6) [1, 2]. The toolkit 
allows simple and fast calculation of power converter circuits including waveform, 
stead-state, and transient analysis using datasheet values of switching devices and 
DC-link capacitors. The aim is to illustrate the influence of PWM control methods 
to students and the interested reader without the use of commercial tools like SIMULINK, 
PLECs or LTSpice. It is clear, that this toolkit cannot be anywhere near the capabilities 
of commercial software but will hopefully provide a better understanding due to the freely 
available source code. The toolkit is obviously not complete; thus, suggestions are 
always welcome.

    
# Publication
The PyPowerSim toolkit is part of the following survey paper and tries to replicate 
the presented architectures and approaches. Please cite the following paper when 
using the PyPowerSim toolkit. When using the B6 architecture and the waveform analysis 
options please also refer to the following article:

Schirmer, Pascal A., Daniel Glose, and Ulrich Ammann. "Zero-voltage and frequency pattern 
selection for DC-link loss minimization in PWM-VSI drives." Electrical Engineering (2022): 1-10.


# Dependencies
The requirements of the PyPowerSim toolkit are summarized in the requirements.txt data file.
In detail, the PyPowerSim Toolkit was implemented using the following dependencies:
- Python 3.8
- Numpy 
- Pandas
- Scipy

The complete list of requirements and packages can be found in "requirements.txt".


# Limitations
Since the toolkit is still under development there are several things that need to be 
improved, are not yet implemented or lack verification with numerical models or measurements.
In the following a list of know issues and limitations is provided:
- The transfer functions for the input and output filter are not yet verified. Also, there is no protection against instability of the transfer functions.
- Soft switching architectures are not included yet.
- The interpolation methods for calculating the tabulated parameter options are only linear now.


# Architecture
The architecture implemented in the PyPowerSim toolkit is exemplary illustrated for a B2 
converter cell in Figure 1. The source code implementation of the PyPowerSim toolkit aims
to follow the data follow of the implementation in Figure 1 for the interested reader to
follow the data and signal flow path through the implementation. The complete description of the
toolkit can be found in \docu.

![img_4.png](docu%2Fimages%2Fimg_4.png)
Figure 1: Proposed converter architecture as implemented in the PyPowerSim toolkit.

As illustrated in Figure 1 the architecture consists of five blocks namely the source, the input
filter, the converter cell, the output filter and the load. Each of these blocks can be freely
configured with the parameter setup described. In the following a short description
without consideration of the filter elements is provided. For the complete description please 
refer to the theory guide and the following theoretical works. 

The converter is power by a voltage source having a constant dc voltage of $V_{dc}$. The converter is controlled by a time
domain switching function $s_a(t)$ translating the constant dc link voltage in a set of high-frequent voltage pulses as
described in (1):

$v_{a0}(t) = s_{a}(t) \cdot \frac{V_{dc}}{2}$     (1)

The switching function can hereby be described as a set of on- and off-states of the high-side switch $S_{a}^{+}$ and the 
low-side switch $S_{a}^{-}$ and can be expressed in the time-domain by the Fourier series in (2):

$s_{a}(t) = \sum_{v=1}^{\infty} \left(a_v cos(v \omega_{el} t) +  b_v cos(v \omega_{el} t)  \right) = \sum_{v=1}^{\infty} c_v e^{-j v \omega_{el} t} $ (2)

where $a_v$,$b_v$,$c_v$ are the Fourier series coefficients, $v$ is the harmonic number and $\omega_{el}$ is the 
electrical circular frequency of the output current. The relation between the current and voltage on the load side can 
then be expressed by the following differential equation in (3):

$v_{L}(t) = R_{L} i_{L}(t) + L_{L} \frac{di_{L}}{dt} + e(t) $ (3)

where $R_{L}$,$L_{L}$ are the resistance and the inductance of the load, $e(t)$ is the induced voltage, and $v_{L}(t)$, $i_{L}(t)$ 
are the load voltage and current respectively.


# Results
In the following chapter a set of reference results is provided using the B6 converter 
architecture and the default setup file. For the default operation an IFX switch is chosen (
IKQ75N120CS6) For a first test run use start.py to calculate the results presented below.

## Model Inputs
In this Section the utilized models are presented. In detail, it includes the transfer functions, the loss models for
semiconductors and capacitors, as well as the reduced order thermal models (transient thermal impedance curves).

![bode.png](docu%2Fimages%2Fbode.png)
Figure 1a: Bode plots for amplitude and phase considering the load, the dc-link, as well as input and output filter.

![model.png](docu%2Fimages%2Fmodel.png)
Figure 1b: Loss models for the semiconductor switches and the capacitor.

![thermal.png](docu%2Fimages%2Fthermal.png)
Figure 1c: Thermal models for the semiconductor switches and the capacitor.

## Sweeping Operation
Below the simulation results of the waveform analysis are displayed. The modulation function is illustrated in Figure 2,
the load currents and voltages in the time-, frequency- and modulation-domain are illustrated in Figure 3 and Figure 4.

![img.png](docu%2Fimages%2Fimg.png)
Figure 2: Modulation function of a B6 converter using space vector modulation and a standard switching sequence 
(0127) with a pulse number of $Q$=21 and $M_i$=1.0

![img_1.png](docu%2Fimages%2Fimg_1.png)
Figure 3: Currents of a B6 converter using space vector modulation and a standard switching sequence (0127) with a 
pulse number of $Q$=21 and $M_i$=1.0. The load angle is equal to $\phi$=17 deg.

![img_2.png](docu%2Fimages%2Fimg_2.png)
Figure 4: Voltages of a B6 converter using space vector modulation and a standard switching sequence (0127) with a 
pulse number of $Q$=21 and $M_i$=1.0. The load angle is equal to $\phi$=17 deg.

## Steady-State Operation
In this Section the results for the steady-state analysis are presented, the results are calculated in closed-loop 
condition, such that the losses are extracted for the stabilized temperature of the junction. The time-domain results
for the currents, the voltages, as well as the conduction and switching losses for the six switches are illustrated 
in Figure 5.

![steady2.png](docu%2Fimages%2Fsteady2.png)
Figure 5: Time-domain currents, voltages, and losses of a B6 converter using space vector modulation and a standard
switching sequence (0127) with a pulse number of $Q$=21 and $M_i$=1.0. The load angle is equal to $\phi$=17 deg.

![steady1.png](docu%2Fimages%2Fsteady1.png)
Figure 6: Time-domain losses and temperatures for the six switches and diodes of a B6 converter using space vector 
modulation and a standard switching sequence (0127) with a pulse number of $Q$=21 and $M_i$=1.0. The load angle is 
equal to $\phi$=17 deg.

## Transient Operation
In this Section the results for the steady-state analysis are presented, the results are calculated in closed-loop 
condition such that the parameters are updated after every fundamental period, i.e. $T_s$=20 ms. The results are averaged
once over the switching sequence, thus displaying junction temperature swing in Figure 7, and once are average over the
fundamental cycle thus illustrating the self-heating due to the internal losses in Figure 8. It should also be noted, 
that in Figure 8 the parameter updates and the positive coupling of the channel resistances and voltages with the 
junction temperature are visible.

![trans2.png](docu%2Fimages%2Ftrans2.png)
Figure 7: Transient time-domain losses and temperatures for the six switches and diodes of a B6 converter using space 
vector modulation and a standard switching sequence (0127) with a pulse number of $Q$=21 and $M_i$=1.0. 
The values are average over the switching cycle.


# Comparison
To get an indication how well the toolkit compares to commercially available solvers several
comparisons with Simulink and PLECs have been conducted (see \docu). The results for comparing to
a Fuji Dual Pack IGBT (2MBI300XBE120) are presented below. The results can be reproduced using the compareFuji.py setup
file under \setup.

| Losses           | PLECs (W) | PyPowerSim (W) | Error (W) | Error (%) | 
|------------------|-----------|----------------|-----------|-----------|
| Transistor (Swi) | 11.1      | 11.5           | 0.40      | 3.60      |
| Transistor (Con) | 18.0      | 17.4           | 0.60      | 3.33      |
| Diode (Swi)      | 5.79      | 6.17           | 0.38      | 6.56      |
| Diode (Con)      | 19.5      | 19.1           | 0.40      | 2.05      |
| Total            | 54.4      | 54.2           | 0.20      | 3.89      |


# Development
As failure and mistakes are inextricably linked to human nature, the toolkit is obviously not perfect, 
thus suggestions and constructive feedback are always welcome. If you want to contribute to the PyPowerSim 
toolkit or spotted any mistake, please contact me via: p.schirmer@herts.ac.uk


# License
The software framework is provided under the MIT license.

# Version History
1) v.0.0: (01.04.2023) Initial version of PyPowerSim
2) v.0.1: (16.06.2023) Major Update and housekeeping
    - Adding MOSFET devices
    - Revising switching loss calculation
    - Adding output control
3) v.0.2: (14.08.2023) Second major review and housekeeping
    - Adding thermal coupling
    - Bugfix switching sequence
    - General housekeeping
4) v.0.3: (13.10.2023) Third major review and housekeeping
    - Fixing steady-state simulation
    - General housekeeping
5) v.0.4: (01.12.2023) Fourth major review and housekeeping
    - Fixing bug in switching losses
    - Adding direct comparison with semiconductor datasheets
    - General housekeeping
   
# References
[1] Holmes, D. Grahame, and Thomas A. Lipo. Pulse width modulation for power converters: principles and practice. 
Vol. 18. John Wiley & Sons, 2003.

[2] Jenni, Felix, and Dieter Wüest. Steuerverfahren für selbstgeführte Stromrichter. vdf Hochschulverlag AG, 1995.
