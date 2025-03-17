#include "Constants.hpp"

/**
 * @namespace Constants
 * @brief Contains definitions of global constants for battery simulations.
 */
namespace Constants {
    // const char* mesh_file = "Mesh_3x90_T3.mesh";            ///< Path to the mesh file
    // const char* dsF_file = "dsF_3x90_T3.txt";               ///< Path to the surface flux data file

    // const char* mesh_file = "Mesh_40x30_A1.mesh";            ///< Path to the mesh file
    // const char* dsF_file = "dsF_40x30_A1.txt";               ///< Path to the surface flux data file

    // const char* mesh_file = "Mesh_3x90_2.mesh";            ///< Path to the mesh file
    // const char* dsF_file = "dsF_3x90_2.txt";               ///< Path to the surface flux data file

    // const char* mesh_file = "Mesh_20x60_2.mesh";            ///< Path to the mesh file
    // const char* dsF_file = "dsF_20x60_2.txt";  

    // const char* mesh_file = "Mesh_40x30_3.mesh";            ///< Path to the mesh file
    // const char* dsF_file = "dsF_40x30_3.txt"; 

    const char* mesh_file = "../II_1_bin.tif";            ///< Path to the mesh file
    const char* dsF_file = "../dsF.txt";

    bool visualization = true;
    const int order = 1;                                    ///< Order of the finite element basis functions
    // const double dh = 0.5e-4;                               ///< Mesh element size // 20x60
    // const double zeta = 1;                        ///< Interfacial thickness // 20x60

    const double dh = 0.2e-4;                               ///< Mesh element size // 3x90
    const double zeta = 1.0 * 0.375;                        ///< Interfacial thickness // 3x90
    // const double ze = 1.0e-5;

    // const double dh = 0.1e-4;                               ///< Mesh element size // 40x30 circ
    // const double zeta = 1.0e-5;                        ///< Interfacial thickness ze = zeta * h

    const double thres = 1.0e-3;                            ///< Threshold value for numerical operations
    const double eps = 1.0e-6;                              ///< Small epsilon value for numerical tolerance
    const double dt = 1.864558472553700e-01 / 40;           ///< Time step size
    const double tm = 0.0;                                  ///< Initial simulation time
    const double t_minus = 7.619047619047619e-01;           ///< Transference number (ratio of cation conductivity to total conductivity)
    const double D0 = 0.00489;                              ///< Base diffusivity
    const double Frd = 96485.3365;                          ///< Faraday constant
    const double Cst1 = 1.6021766e-19 / (1.3806488e-23 * 300.0);   ///< Constant derived from the Nernst equation (e/kT for T = 300 K)
    const double alp = 0.5;                                 ///< Symmetry factor for electrochemical kinetics
    const double rho = 0.0501;                              ///< Lithium site density
    const double Cr = 3.0;                                  ///< C-rate for charging/discharging cycles
    const double Vsr = 0.2;                                 ///< Voltage scanning rate
    const double VCut = 2.7;                                ///< Cut-off voltage
}   
