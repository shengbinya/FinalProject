import numpy as np
import matplotlib.pyplot as plt

# Load the data from the output file
data = np.loadtxt("yb174_optimal_simulation_data.txt") # Changed filename
time = data[:, 0]
position = data[:, 1]
velocity = data[:, 2]

# Create the plots
plt.figure(figsize=(12, 6))

# Plot position vs time
plt.subplot(1, 2, 1)
plt.plot(time, position)
plt.xlabel("Time (s)")
plt.ylabel("Position (m)")
plt.title("Position vs Time (Optimal Parameters)") # Added (Optimal Parameters)
plt.grid(True)

# Plot velocity vs time
plt.subplot(1, 2, 2)
plt.plot(time, velocity)
plt.xlabel("Time (s)")
plt.ylabel("Velocity (m/s)")
plt.title("Velocity vs Time (Optimal Parameters)") # Added (Optimal Parameters)
plt.grid(True)

plt.tight_layout()
plt.show()