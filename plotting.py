import matplotlib.pyplot as plt
import numpy as np

# Load data from the simulation results file
data = np.loadtxt('simulation_results.txt', delimiter=',', skiprows=1)
time = data[:, 2]
x = data[:, 3]
y = data[:, 4]
z = data[:, 5]
vx = data[:, 6]
vy = data[:, 7]
vz = data[:, 8]
delta_values = np.unique(data[:, 0])
power_values = np.unique(data[:, 1])

# Constants (defined as in the C++ code)
Gamma_broad = 2 * np.pi * 29.1e6  # Hz
Is_broad = 59.97e-3  # W/cm^2  <- CHANGED TO mW/cm^2 AND CONVERTED TO W/m^2
Is_broad_SI = Is_broad * 1e4 # W/m^2

# Create plots for each combination of Delta and Power
for delta in delta_values:
    for power in power_values:
        # Filter data for the current Delta and Power
        mask = (data[:, 0] == delta) & (data[:, 1] == power)
        time_filtered = time[mask]
        x_filtered = x[mask]
        y_filtered = y[mask]
        z_filtered = z[mask]
        vx_filtered = vx[mask]
        vy_filtered = vy[mask]
        vz_filtered = vz[mask]

        # Calculate detuning to transition frequency and intensity/saturation intensity
        detuning_over_gamma = delta  # Delta is already normalized to Gamma in C++
        intensity_over_saturation = (2 * power / (np.pi * (1e-2) ** 2)) / Is_broad_SI #using w2D = 1e-2
        print(2 * power / (np.pi * (1e-2) ** 2))
        print(Is_broad_SI)

        # Plot the position and velocity over time for the filtered data
        plt.figure(figsize=(12, 8))
        plt.suptitle(
            f'Detuning/Gamma = {detuning_over_gamma:.2f}, I/Is = {intensity_over_saturation:.2f}')

        plt.subplot(2, 3, 1)
        plt.plot(time_filtered, x_filtered)
        plt.xlabel('Time (s)')
        plt.ylabel('X Position (m)')

        plt.subplot(2, 3, 2)
        plt.plot(time_filtered, y_filtered)
        plt.xlabel('Time (s)')
        plt.ylabel('Y Position (m)')

        plt.subplot(2, 3, 3)
        plt.plot(time_filtered, z_filtered)
        plt.xlabel('Time (s)')
        plt.ylabel('Z Position (m)')

        plt.subplot(2, 3, 4)
        plt.plot(time_filtered, vx_filtered)
        plt.xlabel('Time (s)')
        plt.ylabel('X Velocity (m/s)')

        plt.subplot(2, 3, 5)
        plt.plot(time_filtered, vy_filtered)
        plt.xlabel('Time (s)')
        plt.ylabel('Y Velocity (m/s)')

        plt.subplot(2, 3, 6)
        plt.plot(time_filtered, vz_filtered)
        plt.xlabel('Time (s)')
        plt.ylabel('Z Velocity (m/s)')

        plt.tight_layout()
        plt.show()
