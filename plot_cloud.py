import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Load the data from the output file
data = np.loadtxt("yb174_2d_mot_simulation_data.txt")
initial_positions_x = data[:, 0]
initial_velocities_x = data[:, 1]
initial_positions_y = data[:, 2]
initial_velocities_y = data[:, 3]
final_positions_x = data[:, 4]
final_velocities_x = data[:, 5]
final_positions_y = data[:, 6]
final_velocities_y = data[:, 7]

num_atoms = len(initial_positions_x)

# Constants - these should match the C++ code
hbar = 1.0545718e-34
atomMass = 174.0 * 1.66053906660e-27
gamma = 2 * np.pi * 29e6
lambda_val = 398.9e-9
k = 2 * np.pi / lambda_val
Isat = 59.97
muB = 9.274009994e-24
g_e = 1.0
g_g = 0.0
muPrime = (g_e - g_g) * muB
B0 = 1.0e-3  # Tesla/meter

# Laser parameters (x-direction) - these should be the *optimal* values from C++
S0_x1 = 1.0
delta0_x1 = -2.0 * gamma
S0_x2 = 1.0
delta0_x2 = -2.0 * gamma
# Laser parameters (y-direction)
S0_y1 = 1.0
delta0_y1 = -2.0 * gamma
S0_y2 = 1.0;
delta0_y2 = -2.0 * gamma;

# Time step (should match C++ code)
dt = 1e-6

# Doppler and Zeeman force calculation in 2D
def dopplerAndZeemanForce(x, y, vx, vy):
    # Magnetic field (quadrupole approximation)
    Bx = B0 * x
    By = B0 * y
    Br = np.sqrt(Bx * Bx + By * By)

    # Detunings with Doppler and Zeeman shifts
    delta_x1 = delta0_x1 - vx * k - muPrime * Br / hbar
    delta_x2 = delta0_x2 + vx * k - muPrime * Br / hbar
    delta_y1 = delta0_y1 - vy * k - muPrime * Br / hbar
    delta_y2 = delta0_y2 + vy * k - muPrime * Br / hbar

    # Forces from x-direction lasers
    force_x1 = hbar * k * gamma / 2.0 * S0_x1 / (1.0 + S0_x1 + (delta_x1 / gamma) ** 2)
    force_x2 = -hbar * k * gamma / 2.0 * S0_x2 / (1.0 + S0_x2 + (delta_x2 / gamma) ** 2)

    # Forces from y-direction lasers
    force_y1 = hbar * k * gamma / 2.0 * S0_y1 / (1.0 + S0_y1 + (delta_y1 / gamma) ** 2)
    force_y2 = -hbar * k * gamma / 2.0 * S0_y2 / (1.0 + S0_y2 + (delta_y2 / gamma) ** 2)

    fx = force_x1 + force_x2
    fy = force_y1 + force_y2
    return fx, fy

def dopplerAndZeemanForceX(x, y, vx, vy):
    fx, fy = dopplerAndZeemanForce(x, y, vx, vy)
    return fx

def dopplerAndZeemanForceY(x, y, vx, vy):
    fx, fy = dopplerAndZeemanForce(x, y, vx, vy)
    return fy

# Create the plots
plt.figure(figsize=(18, 6))

# Plot initial and final velocity distributions
plt.subplot(1, 3, 1)
plt.hist(initial_velocities_x, bins=20, alpha=0.5, label='Initial Vx')
plt.hist(final_velocities_x, bins=20, alpha=0.5, label='Final Vx')
plt.hist(initial_velocities_y, bins=20, alpha=0.5, label='Initial Vy')
plt.hist(final_velocities_y, bins=20, alpha=0.5, label='Final Vy')
plt.xlabel("Velocity (m/s)")
plt.ylabel("Number of Atoms")
plt.title("Velocity Distribution")
plt.grid(True)
plt.legend()

# Plot initial position distribution
plt.subplot(1, 3, 2)
plt.scatter(initial_positions_x, initial_positions_y, s=10)
plt.xlabel("Position X (m)")
plt.ylabel("Position Y (m)")
plt.title("Initial Positions")
plt.grid(True)

# Plot final position distribution
plt.subplot(1, 3, 3)
plt.scatter(final_positions_x, final_positions_y, s=10)
plt.xlabel("Position X (m)")
plt.ylabel("Position Y (m)")
plt.title("Final Positions")
plt.grid(True)

plt.tight_layout()
plt.show()

# Create animation of atom positions
fig, ax = plt.subplots()
ax.set_xlim(-0.005, 0.005)
ax.set_ylim(-0.005, 0.005)
scat = ax.scatter(initial_positions_x, initial_positions_y, s=10)
plt.xlabel("Position X (m)")
plt.ylabel("Position Y (m)")
plt.title("Atom Positions Over Time")
plt.grid(True)

def animate(frame):
    global initial_positions_x, initial_velocities_x, initial_positions_y, initial_velocities_y
    for i in range(num_atoms):
        # Euler method (more stable with smaller dt)
        acc_x = dopplerAndZeemanForceX(initial_positions_x[i], initial_positions_y[i], initial_velocities_x[i], initial_velocities_y[i]) / atomMass
        initial_velocities_x[i] += dt * acc_x
        initial_positions_x[i] += dt * initial_velocities_x[i]

        acc_y = dopplerAndZeemanForceY(initial_positions_x[i], initial_positions_y[i], initial_velocities_x[i], initial_velocities_y[i]) / atomMass
        initial_velocities_y[i] += dt * acc_y
        initial_positions_y[i] += dt * initial_velocities_y[i]
    scat.set_offsets(np.c_[initial_positions_x, initial_positions_y])
    return (scat,)

ani = animation.FuncAnimation(fig, animate, frames=100, interval=100, blit=True)
ani.save('atom_motion.gif', writer='pillow', fps=30)
plt.show()