#######################################################################################################################
#######################################################################################################################
# Title:        PWM Distortion Toolkit for Standard Topologies
# Topic:        Power Electronics
# File:         plotGen
# Date:         02.12.2023
# Author:       Dr. Pascal A. Schirmer
# Version:      V.0.4
# Copyright:    Pascal Schirmer
#######################################################################################################################
#######################################################################################################################

#######################################################################################################################
# Import libs
#######################################################################################################################
# ==============================================================================
# Internal
# ==============================================================================

# ==============================================================================
# External
# ==============================================================================
import control
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.pylab as pl
import matplotlib
matplotlib.use('TkAgg')

#######################################################################################################################
# Function
#######################################################################################################################
def plotGenTF(mdl, setupPara, setupTopo):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("INFO: Plotting transfer functions")

    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Parameters
    # ==============================================================================
    fel = setupTopo['fel']
    fsw = setupPara['PWM']['fs']

    # ==============================================================================
    # Variables
    # ==============================================================================
    w = np.linspace(setupTopo['fel']/10, setupPara['PWM']['fs']*10, 1000)

    ###################################################################################################################
    # Pre-Processing
    ###################################################################################################################
    G_load = control.tf(mdl['TF']['Load'].num, mdl['TF']['Load'].den)
    G_inp = control.tf(mdl['TF']['Inp'].num, mdl['TF']['Inp'].den)
    G_out = control.tf(mdl['TF']['Out'].num, mdl['TF']['Out'].den)
    G_dc = control.tf(mdl['TF']['DC'].num, mdl['TF']['DC'].den)

    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    # ==============================================================================
    # General
    # ==============================================================================
    gs = gridspec.GridSpec(2, 4)
    pl.figure()
    txt = "Transfer Functions with fundamental frequency: " + "$f_{el}$=" + str(int(setupTopo['fel'])) + \
          " Hz (blue) and switching frequency: " + "$f_{sw}$=" + str(int(setupPara['PWM']['fs'])) + " Hz (green)"
    pl.suptitle(txt, size=18)
    pl.subplots_adjust(hspace=0.35, wspace=0.20, left=0.075, right=0.925, top=0.90, bottom=0.075)

    # ==============================================================================
    # Input
    # ==============================================================================
    # ------------------------------------------
    # Calc
    # ------------------------------------------
    mag, phase, omega = control.bode(G_inp, w, Hz=True, dB=True, Plot=False)
    wc = np.interp(-np.pi/2, np.flipud(phase), np.flipud(omega))
    Kcu = np.interp(wc, omega, mag)

    # ------------------------------------------
    # Amplitude
    # ------------------------------------------
    pl.subplot(gs[0, 0])
    pl.plot(omega/2/np.pi,  20*np.log10(mag))
    pl.xscale("log")
    pl.ylabel("Magnitude (dB)")
    pl.plot(plt.xlim(), [20*np.log10(Kcu), 20*np.log10(Kcu)], 'r--')
    pl.plot([wc, wc], plt.ylim(), 'r--')
    pl.plot([fel, fel], plt.ylim(), 'b--')
    pl.plot([fsw, fsw], plt.ylim(), 'g--')
    pl.title("Magnitude TF Input with Gain = {0:.3g} dB".format(20*np.log10(Kcu)))
    pl.grid(True, which="both", ls="-")

    # ------------------------------------------
    # Phase
    # ------------------------------------------
    pl.subplot(gs[1, 0])
    pl.plot(omega / 2 / np.pi, phase / np.pi * 180)
    pl.xscale("log")
    pl.xlabel('Frequency (Hz)')
    pl.ylabel("Phase (deg)")
    pl.plot(plt.xlim(), [-180, -180], 'r--')
    pl.plot([wc, wc], plt.ylim(), 'r--')
    pl.plot([fel, fel], plt.ylim(), 'b--')
    pl.plot([fsw, fsw], plt.ylim(), 'g--')
    pl.title("Phase TF Input with Freq = {0:.3g} Hz".format(wc))
    pl.grid(True, which="both", ls="-")

    # ==============================================================================
    # Output
    # ==============================================================================
    # ------------------------------------------
    # Calc
    # ------------------------------------------
    mag, phase, omega = control.bode(G_out, w, Hz=True, dB=True, Plot=False)
    wc = np.interp(-np.pi/2, np.flipud(phase), np.flipud(omega))
    Kcu = np.interp(wc, omega, mag)

    # ------------------------------------------
    # Amplitude
    # ------------------------------------------
    pl.subplot(gs[0, 1])
    pl.plot(omega / 2 / np.pi, 20 * np.log10(mag))
    pl.xscale("log")
    pl.plot(plt.xlim(), [20 * np.log10(Kcu), 20 * np.log10(Kcu)], 'r--')
    pl.plot([wc, wc], plt.ylim(), 'r--')
    pl.plot([fel, fel], plt.ylim(), 'b--')
    pl.plot([fsw, fsw], plt.ylim(), 'g--')
    pl.title("Magnitude TF Output with Gain = {0:.3g} dB".format(20 * np.log10(Kcu)))
    pl.grid(True, which="both", ls="-")

    # ------------------------------------------
    # Phase
    # ------------------------------------------
    pl.subplot(gs[1, 1])
    pl.plot(omega / 2 / np.pi, phase / np.pi * 180)
    pl.xscale("log")
    pl.xlabel('Frequency (Hz)')
    pl.plot(plt.xlim(), [-180, -180], 'r--')
    pl.plot([wc, wc], plt.ylim(), 'r--')
    pl.plot([fel, fel], plt.ylim(), 'b--')
    pl.plot([fsw, fsw], plt.ylim(), 'g--')
    pl.title("Phase TF Output with Freq = {0:.3g} Hz".format(wc))
    pl.grid(True, which="both", ls="-")

    # ==============================================================================
    # DC-Link
    # ==============================================================================
    # ------------------------------------------
    # Calc
    # ------------------------------------------
    mag, phase, omega = control.bode(G_dc, w, Hz=True, dB=True, Plot=False)
    wc = np.interp(-np.pi/2, np.flipud(phase), np.flipud(omega))
    Kcu = np.interp(wc, omega, mag)

    # ------------------------------------------
    # Amplitude
    # ------------------------------------------
    pl.subplot(gs[0, 2])
    pl.plot(omega / 2 / np.pi, 20 * np.log10(mag))
    pl.xscale("log")
    pl.plot(plt.xlim(), [20 * np.log10(Kcu), 20 * np.log10(Kcu)], 'r--')
    pl.plot([wc, wc], plt.ylim(), 'r--')
    pl.plot([fel, fel], plt.ylim(), 'b--')
    pl.plot([fsw, fsw], plt.ylim(), 'g--')
    pl.title("Magnitude TF DC-Link with Gain = {0:.3g} dB".format(20 * np.log10(Kcu)))
    pl.grid(True, which="both", ls="-")

    # ------------------------------------------
    # Phase
    # ------------------------------------------
    pl.subplot(gs[1, 2])
    pl.plot(omega / 2 / np.pi, phase / np.pi * 180)
    pl.xscale("log")
    pl.xlabel('Frequency (Hz)')
    pl.plot(plt.xlim(), [-180, -180], 'r--')
    pl.plot([wc, wc], plt.ylim(), 'r--')
    pl.plot([fel, fel], plt.ylim(), 'b--')
    pl.plot([fsw, fsw], plt.ylim(), 'g--')
    pl.title("Phase TF DC-Link with Freq = {0:.3g} Hz".format(wc))
    pl.grid(True, which="both", ls="-")

    # ==============================================================================
    # Load
    # ==============================================================================
    # ------------------------------------------
    # Calc
    # ------------------------------------------
    mag, phase, omega = control.bode(G_load, w, Hz=True, dB=True, Plot=False)
    wc = np.interp(-np.pi/2, np.flipud(phase), np.flipud(omega))
    Kcu = np.interp(wc, omega, mag)

    # ------------------------------------------
    # Amplitude
    # ------------------------------------------
    pl.subplot(gs[0, 3])
    pl.plot(omega / 2 / np.pi, 20 * np.log10(mag))
    pl.xscale("log")
    pl.plot(plt.xlim(), [20 * np.log10(Kcu), 20 * np.log10(Kcu)], 'r--')
    pl.plot([wc, wc], plt.ylim(), 'r--')
    pl.plot([fel, fel], plt.ylim(), 'b--')
    pl.plot([fsw, fsw], plt.ylim(), 'g--')
    pl.title("Magnitude TF Load with Gain = {0:.3g} dB".format(20 * np.log10(Kcu)))
    pl.grid(True, which="both", ls="-")

    # ------------------------------------------
    # Phase
    # ------------------------------------------
    pl.subplot(gs[1, 3])
    pl.plot(omega / 2 / np.pi, phase / np.pi * 180)
    pl.xscale("log")
    pl.xlabel('Frequency (Hz)')
    pl.plot(plt.xlim(), [-180, -180], 'r--')
    pl.plot([wc, wc], plt.ylim(), 'r--')
    pl.plot([fel, fel], plt.ylim(), 'b--')
    pl.plot([fsw, fsw], plt.ylim(), 'g--')
    pl.title("Phase TF Load with Freq = {0:.3g} Hz".format(wc))
    pl.grid(True, which="both", ls="-")


#######################################################################################################################
# Function
#######################################################################################################################
def plotGenLoss(para, setupPara, setupData):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("INFO: Plotting loss models")

    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    # ==============================================================================
    # Parameters
    # ==============================================================================
    # ------------------------------------------
    # General
    # ------------------------------------------
    Vdc = setupData['stat']['Vdc']
    idx = ['Vf', 'Vfd', 'Eon', 'Eoff', 'Erec', 'Ciss', 'Coss', 'Crss']

    # ------------------------------------------
    # Switch
    # ------------------------------------------
    Vft_con = para['Swi']['Elec']['con']['Vf']
    RonT = para['Swi']['Elec']['con']['Ron']
    Vfd_con = para['Swi']['Elec']['con']['Vfd']
    RonD = para['Swi']['Elec']['con']['RonD']
    Vnom = para['Swi']['Elec']['con']['Vnom']
    Inom = para['Swi']['Elec']['con']['Inom']

    # ==============================================================================
    # Variables
    # ==============================================================================
    # ------------------------------------------
    # Switch
    # ------------------------------------------
    Tj = para['Swi']['Elec']['vec']['Tj'].to_numpy()
    Ift = para['Swi']['Elec']['vec']['If'].to_numpy()
    Ifd = para['Swi']['Elec']['vec']['If'].to_numpy()
    Vf = para['Swi']['Elec']['vec']['Vf'].to_numpy()

    # ------------------------------------------
    # Capacitor
    # ------------------------------------------
    f = para['Cap']['Elec']['vec']['f'].to_numpy()

    ###################################################################################################################
    # Pre-Processing
    ###################################################################################################################
    # ==============================================================================
    # Semiconductor
    # ==============================================================================
    # ------------------------------------------
    # Constant
    # ------------------------------------------
    if setupPara['Elec']['SwiMdl'] == "con":
        # General
        Ciss = para['Swi']['Elec']['con']['Ciss'] * np.ones(np.size(Ift))
        Coss = para['Swi']['Elec']['con']['Coss'] * np.ones(np.size(Ift))
        Crss = para['Swi']['Elec']['con']['Crss'] * np.ones(np.size(Ifd))

        # IGBT
        if setupPara['Elec']['SwiType'] == "IGBT":
            VfT = Vft_con * np.ones(np.size(Ift))
            VfD = Vfd_con * np.ones(np.size(Ifd))
            Eon = para['Swi']['Elec']['con']['Eon'] * np.ones(np.size(Ift))
            Eoff = para['Swi']['Elec']['con']['Eoff'] * np.ones(np.size(Ift))
            Erec = para['Swi']['Elec']['con']['Erec'] * np.ones(np.size(Ifd))

        # MOSFET
        else:
            VfT = Vft_con * np.ones(np.size(Ift))
            VfD = Vfd_con * np.ones(np.size(Ifd))
            Eon = 0.50 * Vnom * Inom * (para['Swi']['Elec']['con']['tr'] + para['Swi']['Elec']['con']['tf']) + 0.5 * para['Swi']['Elec']['con']['Coss'] * Vnom ** 2 * np.ones(np.size(Ift))
            Eoff = 0.50 * Vnom * Inom * (para['Swi']['Elec']['con']['tr'] + para['Swi']['Elec']['con']['tf']) + 0.5 * para['Swi']['Elec']['con']['Coss'] * Vnom ** 2 * np.ones(np.size(Ift))
            Erec = 0.25 * para['Swi']['Elec']['con']['Qrr'] * Vnom * np.ones(np.size(Ifd)) + 0.5 * para['Swi']['Elec']['con']['Crss'] * Vnom ** 2 * np.ones(np.size(Ifd))

    # ------------------------------------------
    # Piece-wise linear (tbi)
    # ------------------------------------------
    elif setupPara['Elec']['SwiMdl'] == "pwl":
        # General
        Ciss = para['Swi']['Elec']['con']['Ciss'] * np.ones(np.size(Ift))
        Coss = para['Swi']['Elec']['con']['Coss'] * np.ones(np.size(Ift))
        Crss = para['Swi']['Elec']['con']['Crss'] * np.ones(np.size(Ifd))

        # IGBT
        if setupPara['Elec']['SwiType'] == "IGBT":
            VfT = Vft_con + RonT * np.abs(Ift)
            VfD = Vfd_con + RonD * np.abs(Ifd)
            Eon = para['Swi']['Elec']['con']['Eon'] * (Ift / Inom) * (Vdc / Vnom)
            Eoff = para['Swi']['Elec']['con']['Eoff'] * (Ift / Inom) * (Vdc / Vnom)
            Erec = para['Swi']['Elec']['con']['Erec'] * (Ifd / Inom) * (Vdc / Vnom)

        # MOSFET
        else:
            VfT = RonT * np.abs(Ift)
            VfD = Vfd_con + RonD * np.abs(Ifd)
            Eon = 0.50 * Vdc * Ift * (para['Swi']['Elec']['con']['tr'] + para['Swi']['Elec']['con']['tf']) + para['Swi']['Elec']['con']['Qrr'] * Vdc
            Eoff = 0.50 * Vdc * Ift * (para['Swi']['Elec']['con']['tr'] + para['Swi']['Elec']['con']['tf'])
            Erec = 0.25 * para['Swi']['Elec']['con']['Qrr'] * Vdc * np.ones(np.size(Ifd))

    # ------------------------------------------
    # Tabular
    # ------------------------------------------
    else:
        VfT = para['Swi']['Elec']['tab']['Vf'].to_numpy()
        VfD = para['Swi']['Elec']['tab']['Vfd'].to_numpy()
        Eon = para['Swi']['Elec']['tab']['Eon'].to_numpy()
        Eoff = para['Swi']['Elec']['tab']['Eoff'].to_numpy()
        Erec = para['Swi']['Elec']['tab']['Erec'].to_numpy()
        Ciss = para['Swi']['Elec']['tab']['Ciss'].to_numpy()
        Coss = para['Swi']['Elec']['tab']['Coss'].to_numpy()
        Crss = para['Swi']['Elec']['tab']['Crss'].to_numpy()

    # ==============================================================================
    # Capacitor
    # ==============================================================================
    # ------------------------------------------
    # Constant
    # ------------------------------------------
    if setupPara['Elec']['CapMdl'] == "con" or setupPara['Elec']['CapMdl'] == "pwl":
        ESR = para['Cap']['Elec']['con']['ESR'] * np.ones(np.size(f))
        tan = para['Cap']['Elec']['con']['tan'] * np.ones(np.size(f))
        C = para['Cap']['Elec']['con']['C'] * np.ones(np.size(f))

    # ------------------------------------------
    # Tabular
    # ------------------------------------------
    else:
        ESR = para['Cap']['Elec']['tab']['ESR']
        tan = para['Cap']['Elec']['tab']['tan']
        C = para['Cap']['Elec']['tab']['C']

    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    # ==============================================================================
    # General
    # ==============================================================================
    gs = gridspec.GridSpec(3, 3)
    pl.figure()
    txt = "Electrical Losses for Switching devices and Capacitors"
    pl.suptitle(txt, size=18)
    pl.subplots_adjust(hspace=0.35, wspace=0.20, left=0.075, right=0.925, top=0.90, bottom=0.075)

    # ==============================================================================
    # Semiconductor
    # ==============================================================================
    # ------------------------------------------
    # Transistor
    # ------------------------------------------
    # Forward Voltage
    try:
        pl.subplot(gs[0, 0])
        pl.plot(Ift, VfT)
        pl.ylabel('$V_{ce/ds}$ (V)')
        pl.xlabel('$I_{ce/ds}$ (A)')
        pl.title("Forward Voltage Transistor")
        pl.grid("on")
    except:
        print("WARN: No tabulated Forward Voltage found")

    # On Switching Energies
    try:
        pl.subplot(gs[1, 0])
        pl.plot(Ift, Eon*1000)
        pl.ylabel('$E_{on}$ (mJ)')
        pl.xlabel('$I_{ce/ds}$ (A)')
        pl.title("On Switching Energies")
        pl.grid("on")
    except:
        print("WARN: No tabulated On Switching Energies found")

    # Off Switching Energies
    try:
        pl.subplot(gs[2, 0])
        pl.plot(Ift, Eoff * 1000)
        pl.ylabel('$E_{off}$ (mJ)')
        pl.xlabel('$I_{ce/ds}$ (A)')
        pl.title("Off Switching Energies")
        pl.grid("on")
    except:
        print("WARN: No tabulated Off Switching Energies found")

    # ------------------------------------------
    # Diode
    # ------------------------------------------
    # Forward Voltage
    try:
        pl.subplot(gs[0, 1])
        pl.plot(Ifd, VfD)
        pl.ylabel('$V_{f}$ (V)')
        pl.xlabel('$I_{f}$ (A)')
        pl.title("Forward Voltage Diode")
        pl.grid("on")
    except:
        print("WARN: No tabulated Forward Voltage found")

    # On Switching Energies
    try:
        pl.subplot(gs[1, 1])
        pl.plot(Ift, Erec * 1000)
        pl.ylabel('$E_{rec}$ (mJ)')
        pl.xlabel('$I_{f}$ (A)')
        pl.title("Reverse Recovery Energies")
        pl.grid("on")
    except:
        print("WARN: No tabulated Reverse Recovery Energies found")

    # Capacitance's
    try:
        pl.subplot(gs[2, 1])
        pl.plot(Vf, Ciss[:, 0] * 1e12)
        pl.plot(Vf, Coss[:, 0] * 1e12)
        pl.plot(Vf, Crss[:, 0] * 1e12)
        pl.ylabel('$C_{iss,oss,rss}$ (pF)')
        pl.xlabel('$V_{f}$ (V)')
        pl.title("Dynamic Capacitance's")
        pl.grid("on")
    except:
        print("WARN: No tabulated capacitance's found")

    # ==============================================================================
    # Capacitor
    # ==============================================================================
    # ------------------------------------------
    # Capacitance
    # ------------------------------------------
    try:
        pl.subplot(gs[0, 2])
        pl.plot(f, C*1000)
        pl.ylabel('$C$ (mF)')
        pl.xlabel('$f$ (Hz)')
        pl.title("Capacitance")
        pl.grid("on")
    except:
        print("WARN: No tabulated capacitance's found")

    # ------------------------------------------
    # ESR
    # ------------------------------------------
    try:
        pl.subplot(gs[1, 2])
        pl.plot(f, ESR*1000)
        pl.ylabel('$ESR$ (mOhm)')
        pl.xlabel('$f$ (Hz)')
        pl.title("Series Resistance")
        pl.grid("on")
    except:
        print("WARN: No tabulated ESR found")

    # ------------------------------------------
    # Loss Angle
    # ------------------------------------------
    try:
        pl.subplot(gs[2, 2])
        pl.plot(f, tan)
        pl.ylabel('$\delta$ (p.u.)')
        pl.xlabel('$f$ (Hz)')
        pl.title("Dissipation Factor")
        pl.grid("on")
    except:
        print("WARN: No tabulated Dissipation Factor found")


#######################################################################################################################
# Function
#######################################################################################################################
def plotGenTher(para):
    ###################################################################################################################
    # MSG IN
    ###################################################################################################################
    print("INFO: Plotting thermal models")

    ###################################################################################################################
    # Initialisation
    ###################################################################################################################
    N = 10000

    ###################################################################################################################
    # Pre-Processing
    ###################################################################################################################
    # ==============================================================================
    # Semiconductor
    # ==============================================================================
    # ------------------------------------------
    # Switch
    # ------------------------------------------
    # Tabulated
    if np.sum(para['Swi']['Ther']['vec']['Rth_JC'].to_numpy()) != 0:
        Rth_JC = para['Swi']['Ther']['vec']['Rth_JC'].to_numpy()
        Cth_JC = para['Swi']['Ther']['vec']['Cth_JC'].to_numpy()

    # Constant
    elif para['Swi']['Ther']['con']['Rth_JC'].to_numpy() != 0:
        Rth_JC = para['Swi']['Ther']['con']['Rth_JC'].to_numpy()
        Cth_JC = para['Swi']['Ther']['con']['Cth_JC'].to_numpy()

    # Default
    else:
        Rth_JC = 0
        Cth_JC = 1

    # ------------------------------------------
    # Diode
    # ------------------------------------------
    # Tabulated
    if np.sum(para['Swi']['Ther']['vec']['Rth_DC'].to_numpy()) != 0:
        Rth_DC = para['Swi']['Ther']['vec']['Rth_DC'].to_numpy()
        Cth_DC = para['Swi']['Ther']['vec']['Cth_DC'].to_numpy()

    # Constant
    elif para['Swi']['Ther']['con']['Rth_DC'].to_numpy() != 0:
        Rth_DC = para['Swi']['Ther']['con']['Rth_DC'].to_numpy()
        Cth_DC = para['Swi']['Ther']['con']['Cth_DC'].to_numpy()

    # Default
    else:
        Rth_DC = 0
        Cth_DC = 1

    # ------------------------------------------
    # Case
    # ------------------------------------------
    # Tabulated
    if np.sum(para['Swi']['Ther']['vec']['Rth_CA'].to_numpy()) != 0:
        Rth_CA = para['Swi']['Ther']['vec']['Rth_CA'].to_numpy()
        Cth_CA = para['Swi']['Ther']['vec']['Cth_CA'].to_numpy()

    # Constant
    elif para['Swi']['Ther']['con']['Rth_CA'].to_numpy() != 0:
        Rth_CA = para['Swi']['Ther']['con']['Rth_CA'].to_numpy()
        Cth_CA = para['Swi']['Ther']['con']['Cth_CA'].to_numpy()

    # Default
    else:
        Rth_CA = 0
        Cth_CA = 1

    # ==============================================================================
    # Capacitor
    # ==============================================================================
    # ------------------------------------------
    # Core
    # ------------------------------------------
    # Tabulated
    if np.sum(para['Cap']['Ther']['vec']['Rth_JC'].to_numpy()) != 0:
        Rth_Cap_JC = para['Cap']['Ther']['vec']['Rth_JC'].to_numpy()
        Cth_Cap_JC = para['Cap']['Ther']['vec']['Cth_JC'].to_numpy()

    # Constant
    elif para['Cap']['Ther']['con']['Rth_JC'].to_numpy() != 0:
        Rth_Cap_JC = para['Cap']['Ther']['con']['Rth_JC'].to_numpy()
        Cth_Cap_JC = para['Cap']['Ther']['con']['Cth_JC'].to_numpy()

    # Default
    else:
        Rth_Cap_JC = 0
        Cth_Cap_JC = 1

    # ------------------------------------------
    # Case
    # ------------------------------------------
    # Tabulated
    if np.sum(para['Cap']['Ther']['vec']['Rth_CA'].to_numpy()) != 0:
        Rth_Cap_CA = para['Cap']['Ther']['vec']['Rth_CA'].to_numpy()
        Cth_Cap_CA = para['Cap']['Ther']['vec']['Cth_CA'].to_numpy()

    # Constant
    elif para['Cap']['Ther']['con']['Rth_JC'].to_numpy() != 0:
        Rth_Cap_CA = para['Cap']['Ther']['con']['Rth_CA'].to_numpy()
        Cth_Cap_CA = para['Cap']['Ther']['con']['Cth_CA'].to_numpy()

    # Default
    else:
        Rth_Cap_CA = 0
        Cth_Cap_CA = 1

    ###################################################################################################################
    # Calculation
    ###################################################################################################################
    # ==============================================================================
    # General
    # ==============================================================================
    gs = gridspec.GridSpec(2, 2)
    pl.figure()
    txt = "Transient Thermal Impedance for Switching Devices and Capacitors"
    pl.suptitle(txt, size=18)
    pl.subplots_adjust(hspace=0.35, wspace=0.20, left=0.075, right=0.925, top=0.90, bottom=0.075)

    # ==============================================================================
    # Switch
    # ==============================================================================
    # ------------------------------------------
    # Device
    # ------------------------------------------
    # Calc
    tau1 = np.abs(Rth_JC * Cth_JC)
    tau2 = np.abs(Rth_DC * Cth_DC)
    t1 = np.linspace(np.min(tau1) / 10, np.max(tau1) * 10, N)
    t2 = np.linspace(np.min(tau2) / 10, np.max(tau2) * 10, N)
    Zth1 = np.zeros((N, len(tau1)))
    Zth2 = np.zeros((N, len(tau2)))
    for i in range(0, len(tau1)):
        Zth1[:, i] = Rth_JC[i] * (1 - np.exp(-t1/tau1[i]))
    Zth1 = np.sum(Zth1, axis=1)
    for i in range(0, len(tau2)):
        Zth2[:, i] = Rth_DC[i] * (1 - np.exp(-t2/tau2[i]))
    Zth2 = np.sum(Zth2, axis=1)

    # Plot
    pl.subplot(gs[0, 0])
    pl.plot(t1, Zth1)
    pl.plot(t2, Zth2)
    pl.xscale("log")
    pl.yscale("log")
    pl.ylabel('$Z_{th}$ (K/W)')
    pl.xlabel('$t$ (sec)')
    pl.title("Thermal Impedance Switch (Junction-Case)")
    pl.legend(['Transistor', 'Diode'])
    pl.grid("on")

    # ------------------------------------------
    # Case
    # ------------------------------------------
    # Calc
    tau1 = np.abs(Rth_CA * Cth_CA)
    t1 = np.linspace(np.min(tau1) / 10, np.max(tau1) * 10, N)
    Zth1 = np.zeros((N, len(tau1)))
    for i in range(0, len(tau1)):
        Zth1[:, i] = Rth_CA[i] * (1 - np.exp(-t1 / tau1[i]))
    Zth1 = np.sum(Zth1, axis=1)

    # Plot
    pl.subplot(gs[1, 0])
    pl.plot(t1, Zth1)
    pl.xscale("log")
    pl.yscale("log")
    pl.ylabel('$Z_{th}$ (K/W)')
    pl.xlabel('$t$ (sec)')
    pl.title("Thermal Impedance Switch (Case-Ambient)")
    pl.grid("on")

    # ==============================================================================
    # Capacitor
    # ==============================================================================
    # ------------------------------------------
    # Device
    # ------------------------------------------
    # Calc
    tau1 = np.abs(Rth_Cap_JC * Cth_Cap_JC)
    t1 = np.linspace(np.min(tau1) / 10, np.max(tau1) * 10, N)
    Zth1 = np.zeros((N, len(tau1)))
    for i in range(0, len(tau1)):
        Zth1[:, i] = Rth_Cap_JC[i] * (1 - np.exp(-t1 / tau1[i]))
    Zth1 = np.sum(Zth1, axis=1)

    # Plot
    pl.subplot(gs[0, 1])
    pl.plot(t1, Zth1)
    pl.xscale("log")
    pl.yscale("log")
    pl.ylabel('$Z_{th}$ (K/W)')
    pl.xlabel('$t$ (sec)')
    pl.title("Thermal Impedance Capacitor (Core-Case)")
    pl.grid("on")

    # ------------------------------------------
    # Case
    # ------------------------------------------
    # Calc
    tau1 = np.abs(Rth_Cap_CA * Cth_Cap_CA)
    t1 = np.linspace(np.min(tau1) / 10, np.max(tau1) * 10, N)
    Zth1 = np.zeros((N, len(tau1)))
    for i in range(0, len(tau1)):
        Zth1[:, i] = Rth_Cap_CA[i] * (1 - np.exp(-t1 / tau1[i]))
    Zth1 = np.sum(Zth1, axis=1)

    # Plot
    pl.subplot(gs[1, 1])
    pl.plot(t1, Zth1)
    pl.xscale("log")
    pl.yscale("log")
    pl.ylabel('$Z_{th}$ (K/W)')
    pl.xlabel('$t$ (sec)')
    pl.title("Thermal Impedance Capacitor (Case-Ambient)")
    pl.grid("on")
