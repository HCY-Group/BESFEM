#include "Constants.hpp"

/**
 * @namespace Constants
 * @brief Contains definitions of global constants for battery simulations.
 */
namespace Constants {

    const char* mesh_file = "../Inputs/Mesh_96x96_P01.mesh";            ///< Path to the mesh file
    const char* dsF_file = "../Inputs/dsF_97x97_P01.txt";               ///< Path to the surface flux data file
    // const char* mesh_file = "../Code_2D/II_1_bin.tif";            ///< Path to the mesh file
    // const char* dsF_file = "../Code_2D/dsF_p.txt";               ///< Path to the surface flux data file

    bool visualization = true;
    const int order = 1;                                    ///< Order of the finite element basis functions
    
    // const double dh = 0.5e-4;                               ///< Mesh element size
    const double dh = 0.0000325;                               ///< Mesh element size

    // const double zeta = 1.0 * 0.375;                        ///< Interfacial thickness
    const double zeta = 1.0;                        ///< Interfacial thickness

    const double thres = 1.0e-3;                            ///< Threshold value for numerical operations
    
    const double eps = 1.0e-6;                              ///< Small epsilon value for numerical tolerance
    // const double eps = 6.75e-10;                              ///< Small epsilon value for numerical tolerance

    
    // const double dt = 1.864558472553700e-01 / 40;           ///< Time step size
    const double dt = 0.0105625;

    const double tm = 0.0;                                  ///< Initial simulation time
    const double t_minus = 7.619047619047619e-01;           ///< Transference number (ratio of cation conductivity to total conductivity)
    const double D0 = 0.00489;                              ///< Base diffusivity
    
    const double Frd = 96485.3365;                          ///< Faraday constant
    const double Cst1 = 1.6021766e-19 / (1.3806488e-23 * 300.0);   ///< Constant derived from the Nernst equation (e/kT for T = 300 K)
    const double alp = 0.5;                                 ///< Symmetry factor for electrochemical kinetics
    
    // const double rho = 0.0501;                              ///< Lithium site density
    const double rho = 0.0312;
    
    // const double Cr = 3.0;                                  ///< C-rate for charging/discharging cycles
    const double Cr = 0.5;

    // const double Vsr = 0.2;                                 ///< Voltage scanning rate
    const double Vsr = 0.009466;                                 ///< Voltage scanning rate

    // const double VCut = 2.7;                                ///< Cut-off voltage
    const double VCut = 2.7;                                ///< Cut-off voltage

}   
