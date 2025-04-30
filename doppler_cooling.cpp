#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip> // Required for setting precision

// Constants (defined as in the Python code)
const double hbar = 1.0545718e-34;
const double kB = 1.380649e-23;
const double muB = 9.274009994e-24;
const double Da = 1.66053906660e-27;
const double g = 9.8;

// Yb 171 parameters
const double m = 171 * Da;
const double Gamma_broad = 2 * M_PI * 29.1e6;
const double Gamma_narrow = 2 * M_PI * 182.4e3;
const double k_broad = (2 * M_PI) / (399e-9);
const double k_narrow = (2 * M_PI) / (555.8e-9);
const double Is_broad = 59.97e1;
const double Is_narrow = 0.139e1;
const double gJ = 1.0;

// Detuning
const double Delta2D = -3;
const double Delta3D = -2.5;
const double DeltaPush = -0.9;

// B-field gradient
const double Bgrad2D = 80.0e2;
const double Bgrad3D = 50.0e2;

// Beam waist
const double wPush = 3.0e-3;
const double w2D = 1.0e-2;
const double w3D = 0.5e-2;

// Beam power
double P_Push = 10e-6; // Initial value, will be swept
double I_Push;
double P_2DMOT = 40e-3; // Initial value, will be swept
double I_2DMOT;
double P_3DMOT = 20e-3; // Initial value, will be swept
double I_3DMOT;

// Saturation parameters
double s0push;
double s02DMOTX1, s02DMOTX2, s02DMOTY1, s02DMOTY2;
double s03DMOTX1, s03DMOTX2, s03DMOTY1, s03DMOTY2, s03DMOTZ1, s03DMOTZ2;

// Capture condition
const double r_coolingRange2D = hbar * abs(Delta2D) * Gamma_broad / (muB * Bgrad2D);
const double l_coolingRange2D = w2D;
const double r_coolingRange3D = hbar * abs(Delta3D) * Gamma_broad / (muB * Bgrad3D);

// Geometry
const double MOT3dXoffset = 0.0;
const double MOT3dYoffset = -0.003;
const double MOT3dZoffset = 0.43;
const double pushbeamMaxZ = 0.1;
const double diaDiffPumpTube = 2e-3;
const double pushTheta = 7.0 * M_PI / 180.0;

// Excited state population
double r22(double s0, double Delta, double Gamma) {
    return (0.5 * s0) / (1 + s0 + pow(2 * Delta / Gamma, 2));
}

// 2D MOT beam functions
double ax2DMOT(double Gamma, double delta, double k, double v, double x, double y, double z) {
    double delta_minus = delta * Gamma - (k * v) - (x * Bgrad2D * muB * gJ / hbar);
    double delta_plus = delta * Gamma + (k * v) + (x * Bgrad2D * muB * gJ / hbar);
    return exp(-2 * y * y / (w2D * w2D)) * exp(-2 * z * z / (w2D * w2D)) * (hbar * k * Gamma / m) *
           (r22(s02DMOTX1, delta_minus, Gamma) - r22(s02DMOTX2, delta_plus, Gamma));
}

double ay2DMOT(double Gamma, double delta, double k, double v, double x, double y, double z) {
    double delta_minus = delta * Gamma - (k * v) - (y * Bgrad2D * muB * gJ / hbar);
    double delta_plus = delta * Gamma + (k * v) + (y * Bgrad2D * muB * gJ / hbar);
    return exp(-2 * x * x / (w2D * w2D)) * exp(-2 * z * z / (w2D * w2D)) * (hbar * k * Gamma / m) *
           (r22(s02DMOTY1, delta_minus, Gamma) - r22(s02DMOTY2, delta_plus, Gamma));
}

double az2DMOT(double Gamma, double delta, double k, double v, double x, double y, double z) {
    return 0.0;
}

double azPush(double Gamma, double delta, double k, double v, double x, double y, double z) {
    if (z < pushbeamMaxZ) {
        return exp(-2 * x * x / (wPush * wPush)) * exp(-2 * y * y / (wPush * wPush)) * (hbar * k * Gamma / m) *
               r22(s0push, delta - (k * v) / Gamma, Gamma);
    } else {
        return 0.0;
    }
}

// 3D MOT beam functions
double ax3DMOT(double Gamma, double delta, double k, double v, double x, double y, double z) {
    double delta_minus = delta * Gamma - (k * v) - ((x - MOT3dXoffset) * Bgrad3D * muB * gJ / hbar);
    double delta_plus = delta * Gamma + (k * v) + ((x - MOT3dXoffset) * Bgrad3D * muB * gJ / hbar);
    return exp(-2 * (y - MOT3dYoffset) * (y - MOT3dYoffset) / (w3D * w3D)) *
           exp(-2 * (z - MOT3dZoffset) * (z - MOT3dZoffset) / (w3D * w3D)) * (hbar * k * Gamma / m) *
           (r22(s03DMOTX1, delta_minus, Gamma) - r22(s03DMOTX2, delta_plus, Gamma));
}

double ay3DMOT(double Gamma, double delta, double k, double v, double x, double y, double z) {
    double delta_minus = delta * Gamma - (k * v) - ((y - MOT3dYoffset) * Bgrad3D * muB * gJ / hbar);
    double delta_plus = delta * Gamma + (k * v) + ((y - MOT3dYoffset) * Bgrad3D * muB * gJ / hbar);
    return exp(-2 * (x - MOT3dXoffset) * (x - MOT3dXoffset) / (w3D * w3D)) *
           exp(-2 * (z - MOT3dZoffset) * (z - MOT3dZoffset) / (w3D * w3D)) * (hbar * k * Gamma / m) *
           (r22(s03DMOTY1, delta_minus, Gamma) - r22(s03DMOTY2, delta_plus, Gamma));
}

double az3DMOT(double Gamma, double delta, double k, double v, double x, double y, double z) {
    double delta_minus = delta * Gamma - (k * v) - ((z - MOT3dZoffset) * Bgrad3D * muB * gJ / hbar);
    double delta_plus = delta * Gamma + (k * v) + ((z - MOT3dZoffset) * Bgrad3D * muB * gJ / hbar);
    return exp(-2 * (x - MOT3dXoffset) * (x - MOT3dXoffset) / (w3D * w3D)) *
           exp(-2 * (y - MOT3dYoffset) * (y - MOT3dYoffset) / (w3D * w3D)) * (hbar * k * Gamma / m) *
           (r22(s03DMOTZ1, delta_minus, Gamma) - r22(s03DMOTZ2, delta_plus, Gamma));
}

// Overall equation of motion
double axMOT(double v, double x, double y, double z) {
    return ax2DMOT(Gamma_broad, Delta2D, k_broad, v, x, y, z) +
           ax3DMOT(Gamma_narrow, Delta3D, k_narrow, v, x - MOT3dXoffset, y - MOT3dYoffset, z - MOT3dZoffset);
}

double ayMOT(double v, double x, double y, double z) {
    return ay2DMOT(Gamma_broad, Delta2D, k_broad, v, x, y, z) +
           ay3DMOT(Gamma_narrow, Delta3D, k_narrow, v, x - MOT3dXoffset, y - MOT3dYoffset, z - MOT3dZoffset) - g +
           azPush(Gamma_narrow, DeltaPush, k_narrow, v, x, y, z) * sin(pushTheta);
}

double azMOT(double v, double x, double y, double z) {
    return az2DMOT(Gamma_broad, Delta2D, k_broad, v, x, y, z) +
           az3DMOT(Gamma_narrow, Delta3D, k_narrow, v, x - MOT3dXoffset, y - MOT3dYoffset, z - MOT3dZoffset) +
           azPush(Gamma_narrow, DeltaPush, k_narrow, v, x, y, z) * cos(pushTheta);
}

// Define the system of ODEs
void mot_ode_system(double t, double y[], double dydt[]) {
    dydt[0] = y[3]; // x' = vx
    dydt[1] = y[4]; // y' = vy
    dydt[2] = y[5]; // z' = vz
    dydt[3] = axMOT(y[3], y[0], y[1], y[2]); // vx' = ax
    dydt[4] = ayMOT(y[4], y[0], y[1], y[2]); // vy' = ay
    dydt[5] = azMOT(y[5], y[0], y[1], y[2]); // vz' = az
}

int main() {
    // Initial conditions
    double initial_position_x = 0.0;
    double initial_position_y = 0.0;
    double initial_position_z = 0.0;
    double initial_velocity_x = 10.0;
    double initial_velocity_y = 0.0;
    double initial_velocity_z = 0.0;

    double initial_conditions[6] = {
        initial_position_x, initial_position_y, initial_position_z,
        initial_velocity_x, initial_velocity_y, initial_velocity_z
    };

    // Time parameters
    double t_start = 0.0;
    double t_end = 0.001; // Reduced time scale
    int num_points = 1000;
    double t_step = (t_end - t_start) / num_points;

    // Create time vector
    std::vector<double> time(num_points);
    for (int i = 0; i < num_points; ++i) {
        time[i] = t_start + i * t_step;
    }

    // Detuning and Power sweep parameters
    std::vector<double> delta_sweep = {-5, -2.5, 0}; // Example detuning values
    std::vector<double> power_sweep = {10e-3, 20e-3, 30e-3, 40e-3}; // Example power values

    // Output file
    std::ofstream outfile("simulation_results.txt");
    if (!outfile.is_open()) {
        std::cerr << "Could not open file for output" << std::endl;
        return 1;
    }
    outfile << "Delta,Power,t,x,y,z,vx,vy,vz\n"; // Added Delta and Power to header

    for (double delta : delta_sweep) {
        for (double power : power_sweep) {
            // Update detuning and power
            //Delta2D = delta; //Uncomment if you want to sweep this
            //Delta3D = delta; //Uncomment if you want to sweep this
            P_Push = power;
            P_2DMOT = power * 4; //scaled
            P_3DMOT = power * 2; //scaled

            // Update intensities
            I_Push = 2 * P_Push / (M_PI * wPush * wPush);
            I_2DMOT = 2 * P_2DMOT / (M_PI * w2D * w2D);
            I_3DMOT = 2 * P_3DMOT / (M_PI * w3D * w3D);

            // Update saturation parameters
            s0push = I_Push / Is_narrow;
            s02DMOTX1 = I_2DMOT / Is_broad;
            s02DMOTX2 = I_2DMOT / Is_broad;
            s02DMOTY1 = I_2DMOT / Is_broad;
            s02DMOTY2 = I_2DMOT / Is_broad;
            s03DMOTX1 = I_3DMOT / Is_broad;
            s03DMOTX2 = I_3DMOT / Is_broad;
            s03DMOTY1 = I_3DMOT / Is_broad;
            s03DMOTY2 = I_3DMOT / Is_broad;
            s03DMOTZ1 = I_3DMOT / Is_broad;
            s03DMOTZ2 = I_3DMOT / Is_broad;

            // Solve ODE system
            std::vector<std::vector<double>> solution(num_points, std::vector<double>(6));
            for (int i = 0; i < 6; ++i) {
                solution[0][i] = initial_conditions[i];
            }

            double y[6];
            double dydt[6];
            for (int i = 0; i < 6; ++i) {
                y[i] = initial_conditions[i];
            }

            for (int i = 0; i < num_points - 1; ++i) {
                mot_ode_system(time[i], y, dydt);
                for (int j = 0; j < 6; ++j) {
                    y[j] = y[j] + dydt[j] * t_step;
                    solution[i + 1][j] = y[j];
                }
            }

            // Output the data to a file
            for (int i = 0; i < num_points; ++i) {
                outfile << delta << "," << power << ","; // Include Delta and Power
                outfile << std::fixed << std::setprecision(6) << time[i] << ","; // More precision
                for (int j = 0; j < 6; ++j) {
                    outfile << std::fixed << std::setprecision(9) << solution[i][j]; // More precision
                    if (j < 5) outfile << ",";
                }
                outfile << "\n";
            }
             std::cout << "Delta = " << delta << ", Power = " << power << " simulation done" << std::endl;
        }
    }

    outfile.close();
    std::cout << "Data written to simulation_results.txt" << std::endl;

    return 0;
}

