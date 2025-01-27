#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

/**
 * @namespace Constants
 * @brief Defines global constants for battery simulations.
 *
 * This namespace contains constants used across the simulation for parameters
 * such as material properties, simulation settings, and physical constants.
 */

namespace Constants {
    
    extern const char* mesh_file;
    extern const char* dsF_file;
    extern const int order;
    extern const double dh;
    extern const double zeta;
    extern const double ze;
    extern const double thres;
    extern const double eps;
    extern const double dt;
    extern const double tm;
    extern const double t_minus;
    extern const double D0;
    extern const double Frd;
    extern const double Cst1;
    extern const double alp;
    extern const double rho;
    extern const double Cr;
    extern const double Vsr;
    extern const double VCut;
}

#endif // CONSTANTS_HPP
