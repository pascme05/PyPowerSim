#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         plotStat_B2
# Date:         01.14.2023
# Author:       Dr. Pascal A. Schirmer
# Version:      V.0.1
# Copyright:    Pascal Schirmer
#######################################################################################################################
#######################################################################################################################

#######################################################################################################################
# Import libs
#######################################################################################################################
# ==============================================================================
# Internal
# ==============================================================================
from src.general.helpFnc import OoM, thd

# ==============================================================================
# External
# ==============================================================================
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import matplotlib.gridspec as gridspec
from scipy.fft import fft

#######################################################################################################################
# Function
#######################################################################################################################
def plotStat_B2(time, freq, setupPara, setupData, setupTopo, setupExp):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("START: Plotting Stationary B2")
    
    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Init
    # ==============================================================================
    # ------------------------------------------
    # Time
    # ------------------------------------------
    timeSw = time['Sw'] 
    timeAc = time['Ac']
    timeDc = time['Dc']
    timeElec = time['Elec']
    timeLoss = time['Loss']
    timeTher = time['Ther']
    
    # ------------------------------------------
    # Frequency
    # ------------------------------------------
    freqSw = freq['Sw']
    freqAc = freq['Ac'] 
    freqDc = freq['Dc']

    # ==============================================================================
    # Parameters
    # ==============================================================================
    fel = setupTopo['fel']
    fs = setupPara['PWM']['fs']
    fsim = setupExp['fsim']
    Q = int(fs/fel)
    R = setupTopo['R']
    L = setupTopo['L']
    Mi = setupData['stat']['Mi']
    Vdc = setupData['stat']['Vdc']
    phiE = setupTopo['phiE']
    down = setupData['stat']['cyc']

    # ==============================================================================
    # Variables
    # ==============================================================================
    t = timeSw['t'].values
    f = fsim * np.linspace(0, 0.5, int(len(t)/2)) / fel

    ###################################################################################################################
    # Pre-Processing
    ###################################################################################################################
    # ==============================================================================
    # General
    # ==============================================================================
    # ------------------------------------------
    # Limits
    # ------------------------------------------
    K = int(np.round((t[-1]-t[0])*fel))
    start = int((len(t) - 1) / K) * (K - 1)
    ende = len(t)
    
    # ------------------------------------------
    # Change time
    # ------------------------------------------
    t = t[start:ende]
    timeSw = timeSw[:][start:ende]
    for c1 in timeAc:
        timeAc[c1] = timeAc[c1][start:ende]
    for c1 in timeDc:
        timeDc[c1] = timeDc[c1][start:ende]
    for c1 in timeElec:
        for c2 in timeElec[c1]:
            timeElec[c1][c2] = timeElec[c1][c2][start:ende]
    for c1 in timeLoss:
        for c2 in timeLoss[c1]:
            timeLoss[c1][c2] = timeLoss[c1][c2][start:ende]
    for c1 in timeTher:
        for c2 in timeTher[c1]:
            timeTher[c1][c2] = timeTher[c1][c2][start:ende]
    
    # ==============================================================================
    # Load angle Total
    # ==============================================================================
    Y = fft(timeAc['v_a'])
    phiV = np.angle(Y)[1]
    Y = fft(timeAc['i_a'])
    phiI = np.angle(Y)[1]
    phi = phiV - phiI + 2*np.pi
    while phi > 2*np.pi:
        phi = phi - 2*np.pi
    
    # ==============================================================================
    # Load angle RL
    # ==============================================================================
    angZ = math.atan2(2*np.pi*fel*L, R)

    # ==============================================================================
    # Distortion
    # ==============================================================================
    I_ph_thd = thd(timeAc['i_a'], t, (2/np.sqrt(2)), 1)
    V_ph_thd = thd(timeAc['v_a'], t, (2/np.sqrt(2)), 1)
    I_dc_thd = thd(timeDc['i_dc'], t, 1, 0)
    V_dc_thd = thd(timeDc['v_dc'], t, 1, 0)
    
    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    # ==============================================================================
    # Switching Function
    # ==============================================================================
    gs = gridspec.GridSpec(3, 2)
    pl.figure()
    txt = "Modulation Functions for: " + "$M_{i}$=" + str(Mi) + "$ ,Q$=" + str(Q) + ", Sampling: " + str(setupPara['PWM']['samp']) + ", Update: " + str(setupPara['PWM']['upd']) + " and Edge Trigger: " + str(setupPara['PWM']['tri'])
    pl.suptitle(txt, size=18)
    pl.subplots_adjust(hspace=0.35, wspace=0.20, left=0.075, right=0.925, top=0.90, bottom=0.075)
    
    # ------------------------------------------
    # Modulation
    # ------------------------------------------
    ax = pl.subplot(gs[0, :])
    pl.plot(t, timeSw['v_ref']/(Vdc/2), t, timeSw['c'])
    pl.title('Carrier and Reference Waveforms')
    pl.xlabel('t in (sec)')
    pl.ylabel('c(t)/r(t) (p.u)')
    pl.grid('on')
    
    # ------------------------------------------
    # Switching Waveform
    # ------------------------------------------
    # Time-Domain
    ax = pl.subplot(gs[1, 0])
    pl.plot(t, timeSw['s'])
    pl.title('Time-domain Switching Function')
    pl.xlabel('t in (sec)')
    pl.ylabel("$s_{a}(t)$ (p.u)")
    pl.grid('on')
    
    # Freq-Domain
    ax = pl.subplot(gs[2, 0])
    pl.stem(f[::down], freqSw['S'][::down])
    pl.xlim(0, 50)
    pl.title('Frequency-domain Switching Function')
    pl.xlabel("$f/f_{1}$ (Hz/Hz)")
    pl.ylabel("$S_{a}(f)$ (p.u)")
    pl.yscale('log')
    pl.ylim(0.0001, )
    pl.grid('on')
    
    # ------------------------------------------
    # Reference Waveform
    # ------------------------------------------
    # Time-Domain
    ax = pl.subplot(gs[1, 1])
    pl.plot(t, timeSw['xs'])
    pl.title('Time-domain Sampled Reference')
    pl.xlabel('t in (sec)')
    pl.ylabel("$x_{a}(t)$ (p.u)")
    pl.grid('on')
    
    # Freq-Domain
    ax = pl.subplot(gs[2, 1])
    pl.stem(f[::down], freqSw['Xs'][::down])
    pl.xlim(0, 50)
    pl.title('Frequency-domain Sampled Reference')
    pl.xlabel("$f/f_{1}$ (Hz/Hz)")
    pl.ylabel("$X_{a}(f)$ (p.u)")
    pl.yscale('log')
    pl.ylim(0.0001, )
    pl.grid('on')
    
    # ==============================================================================
    # Current
    # ==============================================================================
    plt.figure()
    txt = "Currents B2 Bridge for PWM Controll with: " + "$V_{dc}$=" + str(Vdc) + "V, " + "$M_{i}$=" + str(Mi) + "$ ,Q$=" + str(Q) + ", $\phi_{RL}=$" + str(int(math.degrees(angZ))) + "deg" + ", $\phi_{E}=$" + str(int(phiE)) + "deg" + ", $\phi_{VI}=$" + str(int(math.degrees(phi))) + "deg"
    plt.suptitle(txt, size=18)
    plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)
    
    # ------------------------------------------
    # Phase
    # ------------------------------------------
    # Time
    plt.subplot(2,2,1)
    plt.plot(t, timeAc['i_a'])
    plt.ylabel("$i_{a}(t)$ (A)")
    plt.title('Time-domain Currents AC-Side')
    plt.xlabel('time in (sec)')
    plt.grid('on')
    
    # Frequency
    ax = plt.subplot(2,2,2)
    plt.stem(f[::down], freqAc['I_a'][::down])
    plt.ylabel("$I_{a}(f)$ (A)")
    plt.xlim(0, 50)
    plt.title('Frequency-domain Current AC-Side')
    plt.xlabel("$f/f_{1}$ (Hz/Hz)")
    plt.yscale('log')
    plt.ylim(0.1/OoM(max(freqAc['I_a'])), )
    txt = "THD=" + str(round(I_ph_thd*100,2)) + "%"
    plt.text(0.75, 0.90, txt, bbox=dict(facecolor='tab:blue', alpha=0.5), transform=ax.transAxes)
    plt.grid('on')

    # ------------------------------------------
    # DC-Link
    # ------------------------------------------
    # Time
    plt.subplot(2,2,3)
    plt.plot(t, timeDc['i_d_p'], t, np.mean(timeDc['i_d_p'])*np.ones(np.size(timeDc['i_d_p'])), '--')
    plt.ylabel("$i_{dc}(t)$ (A)")
    plt.title('Time-domain Currents DC-Side')
    plt.xlabel('time in (sec)')
    plt.legend(["$i_{dc}^{+}$", "$I_{dc,avg}^{+}$"], loc='upper right')
    plt.grid('on')
    
    # Frequency
    ax = plt.subplot(2,2,4)
    plt.stem(f[::down], freqDc['I_d_p'][::down])
    plt.ylabel("$I_{dc}(f)$ (A)")
    plt.xlim(0, 50)
    plt.title('Frequency-domain Current DC-Side')
    plt.xlabel("$f/f_{1}$ (Hz/Hz)")
    plt.yscale('log')
    plt.ylim(0.1/OoM(max(freqDc['I_d_p'])), )
    txt = "THD=" + str(round(I_dc_thd*100,2)) + "%"
    plt.text(0.75, 0.90, txt, bbox=dict(facecolor='tab:blue', alpha=0.5), transform=ax.transAxes)
    plt.grid('on')
    
    # ==============================================================================
    # Voltage
    # ==============================================================================
    plt.figure()
    txt = "Voltages B2 Bridge for PWM Controll with: " + "$V_{dc}$=" + str(Vdc) + "V, " + "$M_{i}$=" + str(Mi) + ", $\phi_{RL}=$" + str(int(math.degrees(angZ))) + "deg" + ", $\phi_{E}=$" + str(int(phiE)) + "deg" + ", $\phi_{VI}=$" + str(int(math.degrees(phi))) + "deg"
    plt.suptitle(txt, size=18)
    plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)
    
    # ------------------------------------------
    # Phase
    # ------------------------------------------
    # Time
    plt.subplot(2,2,1)
    plt.plot(t, timeAc['v_a0'], t, timeAc['v_a'], t, timeAc['v_L'], t, timeSw['e'])
    plt.ylabel("$v_{a}(t)$ (V)")
    plt.title('Time-domain Voltages AC-Side')
    plt.xlabel('time in (sec)')
    plt.legend(["$v_{a0}(t)$", "$v_{a}(t)$", "$v_{L}(t)$", "$e(t)$"], loc='upper right')
    plt.grid('on')
    
    # Frequency
    ax = plt.subplot(2,2,2)
    plt.stem(f[::down], freqAc['V_a'][::down])
    plt.ylabel("$V_{a}(f)$ (V)")
    plt.xlim(0, 50)
    plt.title('Frequency-domain Voltages AC-Side')
    plt.xlabel("$f/f_{1}$ (Hz/Hz)")
    plt.yscale('log')
    plt.ylim(0.1/OoM(max(freqAc['V_a'])), )
    txt = "THD=" + str(round(V_ph_thd*100,2)) + "%" 
    plt.text(0.75, 0.90, txt, bbox=dict(facecolor='tab:blue', alpha=0.5), transform=ax.transAxes)
    plt.grid('on')

    # ------------------------------------------
    # DC-Link
    # ------------------------------------------
    # Time
    plt.subplot(2,2,3)
    plt.plot(t, timeDc['v_in'], t, timeDc['v_dc'], t, np.mean(timeDc['v_dc'])*np.ones(np.size(timeDc['v_dc'])), '--')
    plt.ylabel("$v_{dc}(t)$ (V)")
    plt.title('Time-domain Voltages DC-Side')
    plt.xlabel('time in (sec)')
    plt.legend(["$v_{in}$", "$v_{dc}$", "$V_{dc,avg}$"], loc='upper right')
    plt.grid('on')
    
    # Frequency
    ax = plt.subplot(2,2,4)
    plt.stem(f[::down], freqDc['V_dc'][::down])
    plt.ylabel("$V_{dc}(f)$ (A)")
    plt.xlim(0, 50)
    plt.title('Frequency-domain Voltages DC-Side')
    plt.xlabel("$f/f_{1}$ (Hz/Hz)")
    plt.yscale('log')
    plt.ylim(0.1/OoM(max(freqDc['V_dc'])), )
    txt = "THD=" + str(round(V_dc_thd*100,2)) + "%"
    plt.text(0.75, 0.90, txt, bbox=dict(facecolor='tab:blue', alpha=0.5), transform=ax.transAxes)
    plt.grid('on')

    # ==============================================================================
    # Time-domain Switches
    # ==============================================================================
    plt.figure()
    txt = "Time domain switches B2 bridge for PWM control with: " + "$V_{dc}$=" + str(Vdc) + "V, " + "$M_{i}$=" + str(Mi) + ", $\phi_{RL}=$" + str(int(math.degrees(angZ))) + "deg" + ", $\phi_{E}=$" + str(int(phiE)) + "deg" + ", $\phi_{VI}=$" + str(int(math.degrees(phi))) + "deg"
    plt.suptitle(txt, size=18)
    plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)
    
    # ------------------------------------------
    # Switch HS
    # ------------------------------------------
    # Current
    ax = plt.subplot(5,2,1)
    plt.plot(t, timeElec['sw']['S1']['i_T'], t, timeElec['sw']['S1']['i_D'])
    plt.ylabel("$i(t)$ (A)")
    plt.title('Time-domain Currents HS Switch')
    ax.set_xticklabels([])
    plt.legend(["$i_{T1}$", "$i_{D1}$"], loc='upper right')
    plt.grid('on')
    
    # Voltage
    ax = plt.subplot(5,2,3)
    plt.plot(t, timeElec['sw']['S1']['v_T'], t, timeElec['sw']['S1']['v_D'])
    plt.ylabel("$v(t)$ (V)")
    plt.title('Time-domain Voltages HS Switch')
    ax.set_xticklabels([])
    plt.legend(["$v_{T1}$", "$v_{D1}$"], loc='upper right')
    plt.grid('on')
    
    # Losses
    ax = plt.subplot(5,2,5)
    plt.plot(t, timeLoss['sw']['S1']['p_T_c'], t, timeLoss['sw']['S1']['p_D_c'])
    plt.ylabel("$p(t)$ (W)")
    plt.title('Time-domain Conduction Losses HS Switch')
    ax.set_xticklabels([])
    plt.legend(["$p_{T1}$", "$p_{D1}$"], loc='upper right')
    plt.grid('on')
    
    # Losses
    ax = plt.subplot(5,2,7)
    plt.plot(t, timeLoss['sw']['S1']['p_T_s'], t, timeLoss['sw']['S1']['p_D_s'])
    plt.ylabel("$p(t)$ (W)")
    plt.title('Time-domain Switching Losses HS Switch')
    ax.set_xticklabels([])
    plt.legend(["$p_{T1}$", "$p_{D1}$"], loc='upper right')
    plt.grid('on')

    # Temperature
    plt.subplot(5,2,9)
    plt.plot(t, timeTher['sw']['T1'], t, timeTher['sw']['D1'])
    plt.ylabel("$\Theta(t)$ (°C)")
    plt.title('Time-domain Thermal HS Switch')
    plt.xlabel('time in (sec)')
    plt.legend(["$\Theta_{T1}$", "$\Theta_{D1}$"], loc='upper right')
    plt.grid('on')

    # ------------------------------------------
    # Switch LS
    # ------------------------------------------
    # Current
    ax = plt.subplot(5,2,2)
    plt.plot(t, timeElec['sw']['S2']['i_T'], t, timeElec['sw']['S2']['i_D'])
    plt.ylabel("$i(t)$ (A)")
    plt.title('Time-domain Currents LS Switch')
    ax.set_xticklabels([])
    plt.legend(["$i_{T1}$", "$i_{D1}$"], loc='upper right')
    plt.grid('on')
    
    # Voltage
    ax = plt.subplot(5,2,4)
    plt.plot(t, timeElec['sw']['S2']['v_T'], t, timeElec['sw']['S2']['v_D'])
    plt.ylabel("$v(t)$ (V)")
    plt.title('Time-domain Voltages LS Switch')
    ax.set_xticklabels([])
    plt.legend(["$v_{T2}$", "$v_{D2}$"], loc='upper right')
    plt.grid('on')
    
    # Losses
    ax = plt.subplot(5,2,6)
    plt.plot(t, timeLoss['sw']['S2']['p_T_c'], t, timeLoss['sw']['S2']['p_D_c'])
    plt.ylabel("$p(t)$ (W)")
    plt.title('Time-domain Conduction Losses LS Switch')
    ax.set_xticklabels([])
    plt.legend(["$p_{T2}$", "$p_{D2}$"], loc='upper right')
    plt.grid('on')

    # Losses
    ax = plt.subplot(5,2,8)
    plt.plot(t, timeLoss['sw']['S2']['p_T_s'], t, timeLoss['sw']['S2']['p_D_s'])
    plt.ylabel("$p(t)$ (W)")
    plt.title('Time-domain Switching Losses LS Switch')
    ax.set_xticklabels([])
    plt.legend(["$p_{T2}$", "$p_{D2}$"], loc='upper right')
    plt.grid('on')
    
    # Temperature
    plt.subplot(5,2,10)
    plt.plot(t, timeTher['sw']['T2'], t, timeTher['sw']['D2'])
    plt.ylabel("$\Theta(t)$ (°C)")
    plt.title('Time-domain Thermal LS Switch')
    plt.xlabel('time in (sec)')
    plt.legend(["$\Theta_{T2}$", "$\Theta_{D2}$"], loc='upper right')
    plt.grid('on')

    # ==============================================================================
    # Time-domain Capacitor
    # ==============================================================================
    plt.figure()
    txt = "Time domain capacitor B2 bridge for PWM control with: " + "$V_{dc}$=" + str(Vdc) + "V, " + "$M_{i}$=" + str(Mi) + ", $\phi_{RL}=$" + str(int(math.degrees(angZ))) + "deg" + ", $\phi_{E}=$" + str(int(phiE)) + "deg" + ", $\phi_{VI}=$" + str(int(math.degrees(phi))) + "deg"
    plt.suptitle(txt, size=18)
    plt.subplots_adjust(hspace=0.35, wspace=0.35, left=0.075, right=0.925, top=0.90, bottom=0.075)
    
    # ------------------------------------------
    # Capacitor
    # ------------------------------------------
    # Current
    plt.subplot(4,1,1)
    plt.plot(t, timeDc['i_c'])
    plt.ylabel("$i(t)$ (A)")
    plt.title('Time-domain Currents DC-Link Capacitor')
    plt.grid('on')
    
    # Voltage
    plt.subplot(4,1,2)
    plt.plot(t, timeDc['v_dc'])
    plt.ylabel("$v(t)$ (V)")
    plt.title('Time-domain Voltages DC-Link Capacitor')
    plt.grid('on')
    
    # Losses
    plt.subplot(4,1,3)
    plt.plot(t, timeLoss['cap']['C1']['p_L'])
    plt.ylabel("$p(t)$ (W)")
    plt.title('Time-domain Losses DC-Link Capacitor')
    plt.grid('on')
    
    # Temperature
    plt.subplot(4,1,4)
    plt.plot(t, timeTher['cap']['C1'])
    plt.ylabel("$\Theta(t)$ (°C)")
    plt.title('Time-domain Thermal DC-Link Capacitor')
    plt.xlabel('time in (sec)')
    plt.grid('on')
    plt.show()

    ###################################################################################################################
    # Post-Processing
    ###################################################################################################################

    ###################################################################################################################
    # MSG OUT
    ###################################################################################################################
    print("END: Plotting Stationary B2")
    
    ###################################################################################################################
    # Return
    ###################################################################################################################