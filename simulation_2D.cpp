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

// Laser parameters (x-direction)
double S0_x1; // Laser 1, positive kx
double delta0_x1;
double S0_x2; // Laser 2, negative kx
double delta0_x2;

// Laser parameters (y-direction)
double S0_y1; // Laser 3, positive ky
double delta0_y1;
double S0_y2; // Laser 4, negative ky
double delta0_y2;

// Doppler force calculation in 2D
void dopplerForce(double vx, double vy, double &fx, double &fy) {
    // Forces from x-direction lasers
    double delta_x1 = delta0_x1 - vx * k;
    double delta_x2 = delta0_x2 + vx * k;
    double force_x1 = hbar * k * gamma / 2.0 * S0_x1 / (1.0 + S0_x1 + pow(delta_x1 / gamma, 2.0));
    double force_x2 = -hbar * k * gamma / 2.0 * S0_x2 / (1.0 + S0_x2 + pow(delta_x2 / gamma, 2.0));

    // Forces from y-direction lasers
    double delta_y1 = delta0_y1 - vy * k;
    double delta_y2 = delta0_y2 + vy * k;
    double force_y1 = hbar * k * gamma / 2.0 * S0_y1 / (1.0 + S0_y1 + pow(delta_y1 / gamma, 2.0));
    double force_y2 = -hbar * k * gamma / 2.0 * S0_y2 / (1.0 + S0_y2 + pow(delta_y2 / gamma, 2.0));

    fx = force_x1 + force_x2;
    fy = force_y1 + force_y2;
}

// Equations of motion for the Yb-174 atom in 2D
struct State {
    double position_x;
    double velocity_x;
    double position_y;
    double velocity_y;
};

// Function to calculate the derivatives of the state variables
State derivatives(const State &state) {
    State dstate;
    double force_x, force_y;
    dopplerForce(state.velocity_x, state.velocity_y, force_x, force_y);
    dstate.position_x = state.velocity_x;
    dstate.velocity_x = force_x / atomMass;
    dstate.position_y = state.velocity_y;
    dstate.velocity_y = force_y / atomMass;
    return dstate;
}

// RK4 integration step
State rk4Step(const State &state, double dt) {
    State k1 = derivatives(state);
    State k2 = derivatives({
        state.position_x + 0.5 * dt * k1.position_x,
        state.velocity_x + 0.5 * dt * k1.velocity_x,
        state.position_y + 0.5 * dt * k1.position_y,
        state.velocity_y + 0.5 * dt * k1.velocity_y});
    State k3 = derivatives({
        state.position_x + 0.5 * dt * k2.position_x,
        state.velocity_x + 0.5 * dt * k2.velocity_x,
        state.position_y + 0.5 * dt * k2.position_y,
        state.velocity_y + 0.5 * dt * k2.velocity_y});
    State k4 = derivatives({
        state.position_x + dt * k3.position_x,
        state.velocity_x + dt * k3.velocity_x,
        state.position_y + dt * k3.position_y,
        state.velocity_y + dt * k3.velocity_y});

    return {
        state.position_x + dt / 6.0 * (k1.position_x + 2.0 * k2.position_x + 2.0 * k3.position_x + k4.position_x),
        state.velocity_x + dt / 6.0 * (k1.velocity_x + 2.0 * k2.velocity_x + 2.0 * k3.velocity_x + k4.velocity_x),
        state.position_y + dt / 6.0 * (k1.position_y + 2.0 * k2.position_y + 2.0 * k3.position_y + k4.position_y),
        state.velocity_y + dt / 6.0 * (k1.velocity_y + 2.0 * k2.velocity_y + 2.0 * k3.velocity_y + k4.velocity_y)};
}

int main() {
    // Initial conditions
    State state = {0.0, 1.0, 0.0, 1.0}; // Initial position (m), initial velocity (m/s) in x and y
    double time = 0.0;
    double dt = 1e-6;
    double simulationTime = 1e-3;

    // Output data file
    std::ofstream outputFile("yb174_2d_optimal_simulation_data.txt");
    if (!outputFile.is_open()) {
        std::cerr << "Error opening output file." << std::endl;
        return 1;
    }
    outputFile << std::fixed << std::setprecision(12);

    // Sweep parameters
    std::vector<double> S0_values = {0.1, 0.5, 1.0, 5.0, 8.0};
    std::vector<double> delta0_values = {-5.0 * gamma, -2.0 * gamma, -0.5 * gamma, 5.0 * gamma, 2.0 * gamma, 0.5 * gamma};

    double minFinalVelocity = 1e9;
    double optimalS0_x1 = 0.0;
    double optimalDelta0_x1 = 0.0;
    double optimalS0_x2 = 0.0;
    double optimalDelta0_x2 = 0.0;
    double optimalS0_y1 = 0.0;
    double optimalDelta0_y1 = 0.0;
    double optimalS0_y2 = 0.0;
    double optimalDelta0_y2 = 0.0;

    // Perform the sweep
    for (double currentS0_x1 : S0_values) {
        for (double currentDelta0_x1 : delta0_values) {
            for (double currentS0_x2 : S0_values) {
                for (double currentDelta0_x2 : delta0_values) {
                    for (double currentS0_y1 : S0_values) {
                        for (double currentDelta0_y1 : delta0_values) {
                            for (double currentS0_y2 : S0_values) {
                                for (double currentDelta0_y2 : delta0_values) {
                                    S0_x1 = currentS0_x1;
                                    delta0_x1 = currentDelta0_x1;
                                    S0_x2 = currentS0_x2;
                                    delta0_x2 = currentDelta0_x2;
                                    S0_y1 = currentS0_y1;
                                    delta0_y1 = currentDelta0_y1;
                                    S0_y2 = currentS0_y2;
                                    delta0_y2 = currentDelta0_y2;
                                    state = {0.0, 1.0, 0.0, 1.0};
                                    time = 0.0;

                                    // Simulate the motion of the atom using the RK4 method
                                    while (time <= simulationTime) {
                                        state = rk4Step(state, dt);
                                        time += dt;
                                    }

                                    // Check for minimum final velocity (magnitude of 2D velocity)
                                    double finalVelocity = sqrt(pow(state.velocity_x, 2.0) + pow(state.velocity_y, 2.0));
                                    if (finalVelocity < minFinalVelocity) {
                                        minFinalVelocity = finalVelocity;
                                        optimalS0_x1 = currentS0_x1;
                                        optimalDelta0_x1 = currentDelta0_x1;
                                        optimalS0_x2 = currentS0_x2;
                                        optimalDelta0_x2 = currentDelta0_x2;
                                        optimalS0_y1 = currentS0_y1;
                                        optimalDelta0_y1 = currentDelta0_y1;
                                        optimalS0_y2 = currentS0_y2;
                                        optimalDelta0_y2 = currentDelta0_y2;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    std::cout << "Optimal S0_x1: " << optimalS0_x1 << ", Optimal Delta0_x1: " << optimalDelta0_x1 / gamma
              << " gamma, Optimal S0_x2: " << optimalS0_x2 << ", Optimal Delta0_x2: " << optimalDelta0_x2 / gamma
              << " gamma, Optimal S0_y1: " << optimalS0_y1 << ", Optimal Delta0_y1: " << optimalDelta0_y1 / gamma
              << " gamma, Optimal S0_y2: " << optimalS0_y2 << ", Optimal Delta0_y2: " << optimalDelta0_y2 / gamma
              << " gamma, Min Final Speed: " << minFinalVelocity << std::endl;

    // Simulate and save data for the optimal case
    S0_x1 = optimalS0_x1;
    delta0_x1 = optimalDelta0_x1;
    S0_x2 = optimalS0_x2;
    delta0_x2 = optimalDelta0_x2;
    S0_y1 = optimalS0_y1;
    delta0_y1 = optimalDelta0_y1;
    S0_y2 = optimalS0_y2;
    delta0_y2 = optimalDelta0_y2;
    state = {0.0, 1.0, 0.0, 1.0};
    time = 0.0;
    while (time <= simulationTime) {
        outputFile << time << " " << state.position_x << " " << state.velocity_x << " " << state.position_y << " " << state.velocity_y << std::endl;
        state = rk4Step(state, dt);
        time += dt;
    }
    outputFile.close();

    std::cout << "Simulation complete. Optimal parameters data written to yb174_2d_optimal_simulation_data.txt" << std::endl;
    return 0;
}