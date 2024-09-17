#include "Constants.hpp"

namespace Constants {
    const char* mesh_file = "Mesh_3x90_T3.mesh";
    const char* dsF_file = "dsF_3x90_T3.txt";
    const int order = 1;
    const double dh = 0.5e-4;
    const double zeta = 1.0 * 0.375;                              // interfacial thickness
    const double thres = 1.0e-3;                                  // AvP musk
    const double eps = 1.0e-6;                                    // var - epsilon
    const double dt = 1.864558472553700e-01 / 40;                 // dt
    const double tm = 0.0;                                        // time
    const double t_minus = 7.619047619047619e-01;                 // transference number
    const double D0 = 0.00489;                                    // base diffusivity 
    const double Frd = 96485.3365;                                // Faraday constant
    const double Cst1 = 1.6021766e-19 / (1.3806488e-23 * 300.0);
    const double alp = 0.5;
    const double rho = 0.0501;                                    // Li site density
    const double Cr = 3.0;                                        // C-rate
    const double Vsr = 0.2;                                       // Voltage scanning rate
    const double Vcut = 2.7;                                      // Cut-off voltage
}   
