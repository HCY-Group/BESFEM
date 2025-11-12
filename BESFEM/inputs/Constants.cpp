#include "Constants.hpp"

/**
 * @namespace Constants
 * @brief Contains definitions of global constants for battery simulations.
 */
namespace Constants {

    const char* mesh_file = "../inputs/disk_Mesh_80x80x6.mesh";            ///< Path to the mesh file
    const char* dsF_file_C = "../inputs/dsF_3x90_r.txt";               ///< Path to the cathode flux data file
    const char* dsF_file_A = "../inputs/dsF_3x90_r.txt";               ///< Path to the anode flux data file

    bool visualization = true;
    const int order = 1;                                    ///< Order of the finite element basis functions
    const double dh = 0.0000325;                            ///< Mesh element size (disk)
    const double gc = 3.3800e-10 *2;			            ///< gradient coefficient
    const double zeta = 1.0;                                ///< Interfacial thickness
    const double thres = 1.0e-3;                            ///< Threshold value for numerical operations
    const double eps = 1.0e-6;                              ///< Small epsilon value for numerical tolerance
    // const double dt = 0.0105625;                            ///< Time step size (HALF CELL)
    const double dt = 0.02 * 3.25e-5 * 3.25e-5 /7.3333e-10;      ///< Time step size (FULL CELL)
    const double tm = 0.0;                                  ///< Initial simulation time
    const double t_minus = 7.619047619047619e-01;           ///< Transference number
    const double D0 = 0.00489;                              ///< Base diffusivity
    const double Frd = 96485.3365;                          ///< Faraday constant
    const double Cst1 = 1.6021766e-19 / (1.3806488e-23 * 300.0);   ///< Constant derived from the Nernst equation (e/kT for T = 300 K)
    const double alp = 0.5;                                 ///< Symmetry factor for electrochemical kinetics
    const double rho_A = 0.0312;                             ///< Anode Lithium site density (graphite)
    const double rho_C = 0.0501;                             ///< Cathode Lithium site density (NMC)
    const double Cr = 2;                                   ///< C-rate for charging/discharging cycles
    const double Vsr0 = 0.009466;                                 ///< Voltage scanning rate (same value for anode, cathode, and electrolyte)
    const double VCut = 0.0;                                ///< Cut-off voltage
    const double init_CnA = 0.95;                        ///< initial concentration in the anode (full)
    const double init_CnC = 0.30;                            ///< initial concentration in the cathode
    const double init_CnE = 0.001;                           ///< initial concentration in the electrolyte
    const double init_BvA = -0.01;                            ///< boundary condition for anode potential (full)
    const double init_BvC = 3.96;                         ///< boundary condition for cathode potential
    const double init_BvE = -0.1;                         ///< boundary condition for electrolyte potential (full & half cathode)

    // constants for half cell - anode side 
    // const double init_CnA = 2.0e-2;                        ///< initial concentration in the anode (half)
    // const double init_BvA = -0.1;                            ///< boundary condition for anode potential (half)
    // const double init_BvE = -0.4686;                               ///< boundary condition for eletrolyte potential (half anode)
    // const double init_CnE = 0.001;                           ///< initial concentration in the electrolyte


    const double init_Rxn = 0.0;                             ///< initial reaction rate
    const double init_RxA = 0.0;                             ///< initial anode reaction
    const double init_RxC = 0.0;                             ///< initial cathode reaction

}   









    // const double Cr = 0.5;          ///< C-rate for charging/discharging cycles
    // const double init_CnA = 2.02e-2;                            ///< initial concentration in the anode
    // const double init_CnC = 0.3;                              ///< initial concentration in the cathode
    // const double init_CnE = 0.001005;                           ///< initial concentration in the electrolyte
    // const double init_BvA = -0.1;                            ///< boundary condition for anode potential
    // const double init_BvC = 2.9395;                         ///< boundary condition for cathode potential
    // const double init_BvE = -0.4686;                         ///< boundary condition for electrolyte potential (anode)
    // const double init_BvE = -1.0;                         ///< boundary condition for electrolyte potential (cathode)