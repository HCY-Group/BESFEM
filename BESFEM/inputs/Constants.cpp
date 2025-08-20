#include "Constants.hpp"

/**
 * @namespace Constants
 * @brief Contains definitions of global constants for battery simulations.
 */
namespace Constants {

    const char* mesh_file = "../inputs/disk_Mesh_80x80x6.mesh";            ///< Path to the mesh file
    const char* dsF_file = "../inputs/disk_dsF_81x81x7.txt";               ///< Path to the surface flux data file
    // const char* mesh_file = "../Code_2D/II_1_bin.tif";            ///< Path to the mesh file
    // const char* dsF_file = "../Code_2D/dsF_p.txt";               ///< Path to the surface flux data file

    // --------------- Constants for Disk Mesh ----------------
    bool visualization = true;
    const int order = 1;                                    ///< Order of the finite element basis functions
    const double dh = 0.0000325;                            ///< Mesh element size (disk)
    const double gc = 3.3800e-10 *2;			            ///< gradient coefficient
    const double zeta = 1.0;                                ///< Interfacial thickness
    const double thres = 1.0e-3;                            ///< Threshold value for numerical operations
    const double eps = 1.0e-6;                              ///< Small epsilon value for numerical tolerance
    const double dt = 0.0105625;                            ///< Time step size
    const double tm = 0.0;                                  ///< Initial simulation time
    const double t_minus = 7.619047619047619e-01;           ///< Transference number (ratio of cation conductivity to total conductivity)
    const double D0 = 0.00489;                              ///< Base diffusivity
    const double Frd = 96485.3365;                          ///< Faraday constant
    const double Cst1 = 1.6021766e-19 / (1.3806488e-23 * 300.0);   ///< Constant derived from the Nernst equation (e/kT for T = 300 K)
    const double alp = 0.5;                                 ///< Symmetry factor for electrochemical kinetics
    const double rho = 0.0312;                             ///< Lithium site density    
    const double Cr = 0.5;          ///< C-rate for charging/discharging cycles
    const double Vsr = 0.009466;                                 ///< Voltage scanning rate
    const double VCut = 0.0;                                ///< Cut-off voltage



    // // --------------- Constants for Rectangle 3x90 Model ----------------
    // const int order = 1;
	// const double dh = 0.2e-4;
	// const double zeta = 1.0 *0.375;				// interfacial thickness
	// const double thres = 1.0e-3;					// AvP musk
	// const double eps = 1.0e-6;					// var-epsilon			
	// const double dt = 1.864558472553700e-01 /40;	// time step
	// const double tm = 0.0;						// time
    // const double t_minus = 7.619047619047619e-01; // transference number
	// const double D0 = 0.00489; 					// base diffusivity	
	// const double Frd = 96485.3365; 				// Faraday constant
	// const double Cst1 = 1.6021766e-19/(1.3806488e-23*300.0);
	// const double alp = 0.5;	
	// const double rho = 0.0501;		 	// Li site density	
	// const double Cr = 3.0;				// C-rate
	// const double Vsr = 0.2;				// voltage scanning rate
	// const double Vcut = 2.7; 				// cut-off voltage
    // const double gc = 3.3800e-10 *2;

}   
