#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <algorithm>

// Constants for Yb-174 atom
const double hbar = 1.0545718e-34;  // Reduced Planck constant (J s)
const double atomMass = 174.0 * 1.66053906660e-27; // Mass of Yb-174 atom (kg)
const double gamma = 2 * M_PI * 29e6;         // Transition linewidth (Hz)
const double lambda = 398.9e-9;      // Transition wavelength (m)
const double k = 2 * M_PI / lambda; // Wave number (m^-1)
const double Isat = 59.97;         // Saturation intensity (mW/cm^2)

// Laser 1 parameters (positive k)
double S0_1;
double delta0_1;

// Laser 2 parameters (negative k)
double S0_2;
double delta0_2;

// Doppler force calculation
double dopplerForce(double velocity) {
    double delta_1 = delta0_1 - velocity * k; // Doppler shifted detuning for laser 1
    double delta_2 = delta0_2 + velocity * k; // Doppler shifted detuning for laser 2
    double force1 = hbar * k * gamma / 2.0 * S0_1 / (1.0 + S0_1 + pow(delta_1 / gamma, 2.0));
    double force2 = -hbar * k * gamma / 2.0 * S0_2 / (1.0 + S0_2 + pow(delta_2 / gamma, 2.0)); // Negative sign for opposite direction
    return force1 + force2;
}

// Equations of motion for the Yb-174 atom
struct State {
    double position; // Position (m)
    double velocity; // Velocity (m/s)
};

// Function to calculate the derivatives of the state variables
State derivatives(const State& state) {
    State dstate;
    double force = dopplerForce(state.velocity);
    dstate.position = state.velocity;
    dstate.velocity = force / atomMass; // Acceleration = Force / Mass
    return dstate;
}

// RK4 integration step
State rk4Step(const State& state, double dt) {
    State k1 = derivatives(state);
    State k2 = derivatives({
        state.position + 0.5 * dt * k1.position,
        state.velocity + 0.5 * dt * k1.velocity
    });
    State k3 = derivatives({
        state.position + 0.5 * dt * k2.position,
        state.velocity + 0.5 * dt * k2.velocity
    });
    State k4 = derivatives({
        state.position + dt * k3.position,
        state.velocity + dt * k3.velocity
    });

    return {
        state.position + dt / 6.0 * (k1.position + 2.0 * k2.position + 2.0 * k3.position + k4.position),
        state.velocity + dt / 6.0 * (k1.velocity + 2.0 * k2.velocity + 2.0 * k3.velocity + k4.velocity)
    };
}

int main() {
    // Initial conditions
    State state = { 0.0, 1.0 }; // Initial position (m), initial velocity (m/s)
    double time = 0.0;             // Initial time (s)
    double dt = 1e-6;            // Time step (s)
    double simulationTime = 4e-3; // Total simulation time (s)

    // Output data file for the optimal case
    std::ofstream outputFile("yb174_optimal_simulation_data.txt");
    if (!outputFile.is_open()) {
        std::cerr << "Error opening output file." << std::endl;
        return 1;
    }
    outputFile << std::fixed << std::setprecision(12);

    // Sweep parameters
    std::vector<double> S0_values = { 0.1, 0.5, 1.0, 5.0, 8.0 };
    std::vector<double> delta0_values = { -5.0*gamma, -4.0*gamma, -3.0*gamma, -2.0*gamma, -0.5*gamma, 5.0*gamma, 4.0*gamma, 3.0*gamma, 2.0*gamma, 0.5*gamma };

    double minFinalVelocity = 1e9;
    double optimalS0_1 = 0.0;
    double optimalDelta0_1 = 0.0;
    double optimalS0_2 = 0.0;
    double optimalDelta0_2 = 0.0;

    // Perform the sweep
    for (double currentS0_1 : S0_values) {
        for (double currentDelta0_1 : delta0_values) {
            for (double currentS0_2 : S0_values) {
                for (double currentDelta0_2 : delta0_values) {
                    S0_1 = currentS0_1;
                    delta0_1 = currentDelta0_1;
                    S0_2 = currentS0_2;
                    delta0_2 = currentDelta0_2;
                    state = { 0.0, 1.0 };
                    time = 0.0;

                    // Simulate the motion of the atom using the RK4 method
                    while (time <= simulationTime) {
                        state = rk4Step(state, dt);
                        time += dt;
                    }

                    // Check for minimum final velocity
                    if (std::abs(state.velocity) < minFinalVelocity) {
                        minFinalVelocity = std::abs(state.velocity);
                        optimalS0_1 = currentS0_1;
                        optimalDelta0_1 = currentDelta0_1;
                        optimalS0_2 = currentS0_2;
                        optimalDelta0_2 = currentDelta0_2;
                    }
                }
            }
        }
    }

    std::cout << "Optimal S0_1: " << optimalS0_1 << ", Optimal Delta0_1: " << optimalDelta0_1/gamma
              << "gamma, Optimal S0_2: " << optimalS0_2 << ", Optimal Delta0_2: " << optimalDelta0_2/gamma
              << "gamma, Min Final Speed: " << minFinalVelocity << std::endl;

    // Simulate and save data for the optimal case
    S0_1 = optimalS0_1;
    delta0_1 = optimalDelta0_1;
    S0_2 = optimalS0_2;
    delta0_2 = optimalDelta0_2;
    state = { 0.0, 1.0 };
    time = 0.0;
    while (time <= simulationTime) {
        outputFile << time << " " << state.position << " " << state.velocity << std::endl;
        state = rk4Step(state, dt);
        time += dt;
    }
    outputFile.close();

    std::cout << "Simulation complete. Optimal parameters data written to yb174_optimal_simulation_data.txt" << std::endl;
    return 0;
}
