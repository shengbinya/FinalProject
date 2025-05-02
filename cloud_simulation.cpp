#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <random>
#include <algorithm>

// Constants for Yb-174 atom
const double hbar = 1.0545718e-34;  // Reduced Planck constant (J s)
const double atomMass = 174.0 * 1.66053906660e-27; // Mass of Yb-174 atom (kg)
const double gamma = 2 * M_PI * 29e6;         // Transition linewidth (Hz)
const double lambda = 398.9e-9;      // Transition wavelength (m)
const double k = 2 * M_PI / lambda; // Wave number (m^-1)
const double Isat = 59.97;         // Saturation intensity (mW/cm^2)
const double muB = 9.274009994e-24; // Bohr magneton (J/T)
const double g_e = 1.0; //example g-factor excited state
const double g_g = 0.0;  //example g-factor ground state
const double muPrime = (g_e - g_g) * muB; // Reduced magnetic moment

// Laser parameters (x-direction)
double S0_x1;
double delta0_x1;
double S0_x2;
double delta0_x2;

// Laser parameters (y-direction)
double S0_y1;
double delta0_y1;
double S0_y2;
double delta0_y2;

// Magnetic field parameters
double B0 = 1.0e-3; // Example magnetic field gradient (T/m)  --  B = B0 * r  , simplified, assuming quadrupole field

// Doppler and Zeeman force calculation in 2D
void dopplerAndZeemanForce(double x, double y, double vx, double vy, double &fx, double &fy) {
    // Magnetic field (quadrupole approximation)
    double Bx = B0 * x;
    double By = B0 * y;
    double Br = sqrt(Bx * Bx + By * By);

    // Detunings with Doppler and Zeeman shifts
    double delta_x1 = delta0_x1 - vx * k - muPrime * Br / hbar; //for simplicity, assume linear polarization, so Zeeman shift is the same for both lasers
    double delta_x2 = delta0_x2 + vx * k - muPrime * Br / hbar;
    double delta_y1 = delta0_y1 - vy * k - muPrime * Br / hbar;
    double delta_y2 = delta0_y2 + vy * k - muPrime * Br / hbar;

    // Forces from x-direction lasers
    double force_x1 = hbar * k * gamma / 2.0 * S0_x1 / (1.0 + S0_x1 + pow(delta_x1 / gamma, 2.0));
    double force_x2 = -hbar * k * gamma / 2.0 * S0_x2 / (1.0 + S0_x2 + pow(delta_x2 / gamma, 2.0));

    // Forces from y-direction lasers
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
    dopplerAndZeemanForce(state.position_x, state.position_y, state.velocity_x, state.velocity_y, force_x, force_y);
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
    // Simulation parameters
    int numAtoms = 1000;
    double time = 0.0;
    double dt = 1e-6;
    double simulationTime = 1e-2;

    // Output data file
    std::ofstream outputFile("yb174_2d_mot_simulation_data.txt");
    if (!outputFile.is_open()) {
        std::cerr << "Error opening output file." << std::endl;
        return 1;
    }
    outputFile << std::fixed << std::setprecision(12);

    // Random number generation for initial conditions
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> posDistX(0.0, 1e-3); // Mean 0, std dev 1mm
    std::normal_distribution<> posDistY(0.0, 1e-3);
    std::normal_distribution<> velDistX(1.0, 2.0); // Mean 0, std dev 1m/s
    std::normal_distribution<> velDistY(1.0, 2.0);

    // Store initial and final positions and velocities for all atoms
    std::vector<std::vector<double>> initialPositionsX(numAtoms);
    std::vector<std::vector<double>> initialVelocitiesX(numAtoms);
    std::vector<std::vector<double>> initialPositionsY(numAtoms);
    std::vector<std::vector<double>> initialVelocitiesY(numAtoms);
    std::vector<std::vector<double>> finalPositionsX(numAtoms);
    std::vector<std::vector<double>> finalVelocitiesX(numAtoms);
    std::vector<std::vector<double>> finalPositionsY(numAtoms);
    std::vector<std::vector<double>> finalVelocitiesY(numAtoms);

    // Laser parameters (x-direction) - these should be the *optimal* values found from a parameter sweep
    S0_x1 = 2.5;
    delta0_x1 = -1.0 * gamma;
    S0_x2 = 2.5;
    delta0_x2 = -1.0 * gamma;

    // Laser parameters (y-direction)
    S0_y1 = 2.5;
    delta0_y1 = -1.0 * gamma;
    S0_y2 = 2.5;
    delta0_y2 = -1.0 * gamma;

    // Simulate the motion of each atom
    std::vector<State> atoms(numAtoms);
    for (int i = 0; i < numAtoms; ++i) {
        // Initialize atom positions and velocities
        atoms[i].position_x = posDistX(gen);
        atoms[i].velocity_x = velDistX(gen);
        atoms[i].position_y = posDistY(gen);
        atoms[i].velocity_y = velDistY(gen);

        // Store initial conditions
        initialPositionsX[i].push_back(atoms[i].position_x);
        initialVelocitiesX[i].push_back(atoms[i].velocity_x);
        initialPositionsY[i].push_back(atoms[i].position_y);
        initialVelocitiesY[i].push_back(atoms[i].velocity_y);

        // Simulate the motion of the atom using the RK4 method
        time = 0.0;
        while (time <= simulationTime) {
            atoms[i] = rk4Step(atoms[i], dt);
            time += dt;
        }

        // Store final conditions
        finalPositionsX[i].push_back(atoms[i].position_x);
        finalVelocitiesX[i].push_back(atoms[i].velocity_x);
        finalPositionsY[i].push_back(atoms[i].position_y);
        finalVelocitiesY[i].push_back(atoms[i].velocity_y);

        // Write data to file
        outputFile << initialPositionsX[i][0] << " " << initialVelocitiesX[i][0] << " "
                   << initialPositionsY[i][0] << " " << initialVelocitiesY[i][0] << " "
                   << finalPositionsX[i][0] << " " << finalVelocitiesX[i][0] << " "
                   << finalPositionsY[i][0] << " " << finalVelocitiesY[i][0] << std::endl;
    }

    // Close the output file
    outputFile.close();

    std::cout << "Simulation complete. Data written to yb174_2d_mot_simulation_data.txt" << std::endl;
    return 0;
}