#include "../include/Constants.hpp"

/**
 * @namespace Constants
 * @brief Contains definitions of global constants for battery simulations.
 */
namespace Constants {

    const double Frd = 96485.3365;
    const double RT = 2494.33859;
    const double alp = 0.5;
    const double Cst1 = 1.6021766e-19 / (1.3806488e-23 * 300.0);
    const int order = 1;

    const double t_minus = 7.619047619047619e-01;
    const double D0 = 0.00489;

    const double zeta = 1.0;
    const double thres = 1.0e-3;
    const double eps = 1.0e-6;

    const double Perm = 1e-7;

    const double init_Rxn = 1e-7;
    const double init_RxA = 1e-7;
    const double init_RxC = 1e-7;

}

    // bool visualization = true;
    // const int order = 1;                                    ///< Order of the finite element basis functions
    // // const double dh = 2.01005e-05;                            ///< Mesh element size (TIFF)
    // const double dh = 5.0e-06;
    // const double gc = 3.3800e-10 * 3.0 ;			            ///< gradient coefficient
    // const double zeta = 1.0;                                ///< Interfacial thickness
    // const double thres = 1.0e-3;                            ///< Threshold value for numerical operations
    // const double eps = 1.0e-6;                              ///< Small epsilon value for numerical tolerance
    // // const double dt = 0.0105625;                            ///< Time step size 
    // const double dt = 0.001;
    // const double tm = 0.0;                                  ///< Initial simulation time
    // const double t_minus = 7.619047619047619e-01;           ///< Transference number
    // const double D0 = 0.00489;                              ///< Base diffusivity
    // const double Frd = 96485.3365;                          ///< Faraday constant
    // const double Cst1 = 1.6021766e-19 / (1.3806488e-23 * 300.0);   ///< Constant derived from the Nernst equation (e/kT for T = 300 K)
    // const double alp = 0.5;                                 ///< Symmetry factor for electrochemical kinetics
    // const double rho_A = 0.0312;                             ///< Anode Lithium site density (graphite)
    // // const double rho_C = 0.0501;                             ///< Cathode Lithium site density (NMC)
    // const double rho_LFP = 0.02273544498;
    // const double Cr = 1.0;                                   ///< C-rate for charging/discharging cycles
    // // const double Vsr0 = 0.009466;                                 ///< Voltage scanning rate (same value for anode, cathode, and electrolyte)
    // const double VCut = 0.0;                                ///< Cut-off voltage
    // const double Vsr0 = 2.0;

    // // const double init_CnA = 0.95;                        ///< initial concentration in the anode (full)
    // // const double init_CnC = 0.30;                            ///< initial concentration in the cathode
    // // const double init_CnE = 0.001;                           ///< initial concentration in the electrolyte
    // // const double init_BvA = -0.01;                            ///< boundary condition for anode potential (full)
    // // const double init_BvC = 3.96;                         ///< boundary condition for cathode potential
    // // const double init_BvE = -0.1;                         ///< boundary condition for electrolyte potential (full & half cathode)

    // const double rho_C = 0.02273544498;                      ///< Cathode Lithium site density (LFP)
    // const double RT = 2494.33859;                           ///< RT constant at 300K
    // const double Perm = 1e-7;                              ///< Permittivity constant (paper uses 1.0e-7)

    // // // constants for half cell - anode side 
    // // const double init_CnA = 2.0e-2;                        ///< initial concentration in the anode (half)
    // // const double init_BvA = -0.1;                            ///< boundary condition for anode potential (half)
    // // const double init_BvE = -0.4686;                               ///< boundary condition for eletrolyte potential (half anode)
    // // const double init_CnE = 0.001;                           ///< initial concentration in the electrolyte

    // const double init_Rxn = 1e-7;                             ///< initial reaction rate
    // const double init_RxA = 1e-7;                             ///< initial anode reaction
    // const double init_RxC = 1e-7;                             ///< initial cathode reaction

// }   









    // const double Cr = 0.5;          ///< C-rate for charging/discharging cycles
    // const double init_CnA = 2.02e-2;                            ///< initial concentration in the anode
    // const double init_CnC = 0.3;                              ///< initial concentration in the cathode
    // const double init_CnE = 0.001005;                           ///< initial concentration in the electrolyte
    // const double init_BvA = -0.1;                            ///< boundary condition for anode potential
    // const double init_BvC = 2.9395;                         ///< boundary condition for cathode potential
    // const double init_BvE = -0.4686;                         ///< boundary condition for electrolyte potential (anode)
    // const double init_BvE = -1.0;                         ///< boundary condition for electrolyte potential (cathode)