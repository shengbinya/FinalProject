import numpy as np
import matplotlib.pyplot as plt

# Load the data from the output file
#data = np.loadtxt("yb174_2d_optimal_simulation_data.txt")
data = np.loadtxt("yb174_2d_mot_magnetic_simulation_data.txt")
time = data[:, 0]
position_x = data[:, 1]
velocity_x = data[:, 2]
position_y = data[:, 3]
velocity_y = data[:, 4]

# Create the plots
plt.figure(figsize=(12, 6))

# Plot position vs time in x and y
plt.subplot(1, 2, 1)
plt.plot(time, position_x, label='Position X')
plt.plot(time, position_y, label='Position Y')
plt.xlabel("Time (s)")
plt.ylabel("Position (m)")
plt.title("Position vs Time")
plt.grid(True)
plt.legend()

# Plot velocity vs time in x and y
plt.subplot(1, 2, 2)
plt.plot(time, velocity_x, label='Velocity X')
plt.plot(time, velocity_y, label='Velocity Y')
plt.xlabel("Time (s)")
plt.ylabel("Velocity (m/s)")
plt.title("Velocity vs Time")
plt.grid(True)
plt.legend()

"""
# Plot position x vs position y
plt.subplot(2, 2, 3)
plt.plot(position_x, position_y)
plt.xlabel("Position X (m)")
plt.ylabel("Position Y (m)")
plt.title("Position X vs Position Y")
plt.grid(True)

# Plot velocity x vs velocity y
plt.subplot(2, 2, 4)
plt.plot(velocity_x, velocity_y)
plt.xlabel("Velocity X (m/s)")
plt.ylabel("Velocity Y (m/s)")
plt.title("Velocity X vs Velocity Y")
plt.grid(True)
"""

plt.tight_layout()
plt.show()
