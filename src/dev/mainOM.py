import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

###################################################################################################################
# Initialisation
###################################################################################################################
# Define the modulation index and the desired angle (alpha)
Ns = 501
fs = 1/Ns
Mi_t = np.linspace(0, 1.25, 100)  # Modulation index (0 < Mi <= 1)
alpha_t = np.linspace(0, 2 * np.pi, Ns)  # Desired angle sweeping from 0 to 2π

# Define transformation matrices for each sector
Tc = np.zeros((2, 2, 6))
Tc[:, :, 0] = np.array([[1, 0.5], [0, np.sqrt(3) / 2]])
Tc[:, :, 1] = np.array([[0.5, -0.5], [np.sqrt(3) / 2, np.sqrt(3) / 2]])
Tc[:, :, 2] = np.array([[-0.5, -1], [np.sqrt(3) / 2, 0]])
Tc[:, :, 3] = np.array([[-1, -0.5], [0, -np.sqrt(3) / 2]])
Tc[:, :, 4] = np.array([[-0.5, 0.5], [-np.sqrt(3) / 2, -np.sqrt(3) / 2]])
Tc[:, :, 5] = np.array([[0.5, 1], [-np.sqrt(3) / 2, 0]])

# Initialize arrays for alpha and beta voltages
v_m = np.zeros(len(Mi_t))

###################################################################################################################
# Perform SVM
###################################################################################################################
for jj in range(0, len(Mi_t)):
    Mi = Mi_t[jj]
    v_al = np.zeros_like(alpha_t)
    v_be = np.zeros_like(alpha_t)
    v_st = np.zeros_like(alpha_t)
    v_0 = np.zeros_like(alpha_t)
    for j in range(len(alpha_t)):
        # Determine the sector (R) based on the angle alpha
        alpha = alpha_t[j]
        R = int(np.floor(alpha / (np.pi / 3)))

        # ==============================================================================
        # Under-Modulation
        # ==============================================================================
        # Calculate the angle relative to the start of the sector
        alpha = alpha - R * (np.pi / 3)

        # Calculate t1 and t2 based on the sector
        t1 = np.sqrt(3) / 2 * Mi * np.sin(np.pi / 3 - alpha)
        t2 = np.sqrt(3) / 2 * Mi * np.sin(alpha)

        # Calculate the total zero vector time
        Tz = 1 - t1 - t2
        t0 = Tz * 0.5
        t7 = Tz * 0.5

        # ==============================================================================
        # Over-Modulation I
        # ==============================================================================
        if 2 / np.sqrt(3) < Mi <= 0.9517 / np.pi * 4:
            # Angle Intersection
            if 0.9068 / np.pi * 4 <= Mi < 0.9095 / np.pi * 4:
                alpha_r = -30.23 * Mi / 4 * np.pi + 27.94
            elif 0.9095 / np.pi * 4 < Mi < 0.9485 / np.pi * 4:
                alpha_r = -8.58 * Mi / 4 * np.pi + 8.23
            else:
                alpha_r = -26.43 * Mi / 4 * np.pi + 25.15

            # Amplitude Correction
            Mi_om = 2 / (np.sqrt(3) * np.cos(np.pi / 6 - alpha_r))

            # Switching Times
            t1 = np.sqrt(3) / 2 * Mi_om * np.sin(np.pi / 3 - alpha)
            t2 = np.sqrt(3) / 2 * Mi_om * np.sin(alpha)

            # Calculate the total zero vector time
            Tz = 1 - t1 - t2
            t0 = Tz * 0.5
            t7 = Tz * 0.5

            # Negative Zero Vectors
            if Tz < 0:
                t1 = (np.sqrt(3) * np.cos(alpha) - np.sin(alpha)) / (np.sqrt(3) * np.cos(alpha) + np.sin(alpha))
                t2 = 1 - t1
                t0 = 0
                t7 = 0

        # ==============================================================================
        # Over-Modulation II
        # ==============================================================================
        if 0.9517 / np.pi * 4 < Mi < 4 / np.pi:
            # Holding angle
            if 0.9517 / np.pi * 4 <= Mi < 0.9800 / np.pi * 4:
                alpha_h = 6.40 * Mi / 4 * np.pi - 6.09
            elif 0.9800 / np.pi * 4 < Mi < 0.9975 / np.pi * 4:
                alpha_h = 11.75 * Mi / 4 * np.pi - 11.34
            else:
                alpha_h = 48.96 * Mi / 4 * np.pi - 48.43

            # Comparison
            if 0 <= alpha < alpha_h:
                alpha_o = 0
                t1 = 1
                t2 = 0
            elif alpha_h <= alpha < np.pi / 3 - alpha_h:
                alpha_o = (np.pi / 6) * (alpha - alpha_h) / (np.pi / 6 - alpha_h)
                t1 = (np.sqrt(3) * np.cos(alpha_o) - np.sin(alpha_o)) / (np.sqrt(3) * np.cos(alpha_o) + np.sin(alpha_o))
                t2 = 1 - t1
            else:
                alpha_o = np.pi / 3
                t1 = 0
                t2 = 1

            # Times
            t0 = 0
            t7 = 0

        # ==============================================================================
        # Over-Modulation III
        # ==============================================================================
        if Mi == 4 / np.pi:
            if alpha < np.pi / 6:
                t1 = 1
                t2 = 0
            else:
                t2 = 1
                t1 = 0
            t0 = 0
            t7 = 0

        # ==============================================================================
        # Post-Processing
        # ==============================================================================
        # Flip alpha for sectors 0-1, 2-3, 4-5 (which correspond to angles 0-60, 120-180, 240-300 degrees)
        if R in [0, 2, 4]:
            dt_v = np.array([t1, t2])
        else:
            dt_v = np.array([t2, t1])

        # Adjust sector index R for matrix selection
        R = R % 6  # Ensure R is within 0 to 5

        # Calculate the voltage vector in alpha-beta plane
        dt_v = np.array([t1, t2])
        v = np.squeeze(Tc[:, :, R]) @ dt_v
        v_al[j] = 2 / 3 * v[0]
        v_be[j] = 2 / 3 * v[1]
        v_0[j] = 0.5 * (-t0 - (-1) ** (R + 1) * t1 / 3 + (-1) ** (R + 1) * t2 / 3 + t7)
        v_st[j] = np.sqrt(v_al[j] ** 2 + v_be[j] ** 2)

        # FFT of the signal
        N = len(v_al)
        v_a_fft = np.fft.fft(v_al)
        v_a_fft = 2.0 / N * np.abs(v_a_fft[:N // 2])  # Take the magnitude and normalize
        freqs = np.fft.fftfreq(N, 1 / fs)[:N // 2]  # Frequency bins

        # Find the fundamental frequency component
        fundamental_idx = np.argmax(v_a_fft[1:]) + 1  # Skip the zero frequency (DC)
        fundamental_freq = freqs[fundamental_idx]
        v_m[jj] = v_a_fft[fundamental_idx]

###################################################################################################################
# Post-processing
###################################################################################################################
df = pd.DataFrame(v_al, columns=['v_a0'])
df.to_excel('output.xlsx', index=False)

###################################################################################################################
# Plotting
###################################################################################################################
# Plot alpha and beta voltages
plt.plot(Mi_t, v_m)
plt.title("Fundamental of v_a")
plt.xlabel("Mi (p.u.)")
plt.ylabel("Amplitude")
plt.grid(True)

# Plot alpha and beta voltages
plt.figure(figsize=(8, 6))
plt.plot(alpha_t, v_al, label=r'$\alpha$ Voltage')
plt.plot(alpha_t, v_be, label=r'$\beta$ Voltage')
plt.plot(alpha_t, v_0, label=r'Zero Voltage')
plt.plot(alpha_t, v_st, label=r'Ref Voltage')
plt.title('Alpha, Beta and Zero Voltages using SVM')
plt.xlabel('Angle (radians)')
plt.ylabel('Voltage (normalized)')
plt.legend()
plt.grid(True)

# Scaling factor
peak_voltage = 2 / 3  # Given peak value of the voltage vector

# Hexagon vertices (corresponding to active vectors)
hexagon_angles = np.linspace(0, 2 * np.pi, 7)  # 7 points to close the hexagon
hexagon_radius = peak_voltage
hexagon_x = hexagon_radius * np.cos(hexagon_angles)
hexagon_y = hexagon_radius * np.sin(hexagon_angles)

# Inner circle parameters
circle_radius = peak_voltage * (np.sqrt(3) / 2)  # Circle inside hexagon
circle_theta = np.linspace(0, 2 * np.pi, 100)
circle_x = circle_radius * np.cos(circle_theta)
circle_y = circle_radius * np.sin(circle_theta)

# Plotting the space vector hexagon and the inner circle
plt.figure(figsize=(6, 6))
plt.plot(hexagon_x, hexagon_y, 'b-', label='Space Vector Hexagon')
plt.plot(circle_x, circle_y, 'r--', label='Inner Circle')

# Adding labels and formatting the plot
plt.title('Space Vector Hexagon and Inner Circle')
plt.xlabel(r'$\alpha$ (Real Axis)')
plt.ylabel(r'$\beta$ (Imaginary Axis)')
plt.axhline(0, color='black', linewidth=0.5)
plt.axvline(0, color='black', linewidth=0.5)
plt.legend()
plt.grid(True)
plt.axis('equal')  # Ensures equal scaling on both axes

# Plot alpha and beta voltages
plt.plot(v_al, v_be, 'o', label=r'Ref Voltage')
plt.title('Alpha and Beta Voltages using SVM')
plt.xlabel('V_al')
plt.ylabel('V_be')
plt.legend()
plt.grid(True)
plt.show()
