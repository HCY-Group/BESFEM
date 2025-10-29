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
    extern const char* dsF_file_A;
    extern const char* dsF_file_C;
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
    extern const double rho_A;
    extern const double rho_C;
    extern const double Cr;
    extern const double Vsr0;
    extern const double VCut;
    extern const double gc;
    extern const double init_CnA;
    extern const double init_CnC;
    extern const double init_CnE;
    extern const double init_BvA;
    extern const double init_BvC;
    extern const double init_BvE;
    extern const double init_Rxn;
    extern const double init_RxA;
    extern const double init_RxC;
}

#endif // CONSTANTS_HPP
