#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

/**
 * @file Constants.hpp
 * @brief Declares global constants used throughout BESFEM simulations.
 *
 * Includes mesh/input file names, polynomial order, grid spacing parameters,
 * numerical tolerances, electrochemical constants, and initial conditions.
 */

/**
 * @namespace Constants
 * @brief Global user-adjustable simulation constants.
 *
 * This namespace provides a central location for:
 * - mesh and distance-function file paths,
 * - finite element order,
 * - physical parameters,
 * - model constants,
 * - default initialization values for potentials, concentrations, and reactions.
 *
 * Many of these values are read at runtime or overridden by command-line
 * arguments (see `SimulationConfig`).
 */
namespace Constants {

    // -------------------------------------------------------------------------
    // File paths
    // -------------------------------------------------------------------------

    extern const char* mesh_file;   ///< Default mesh file path.
    extern const char* dsF_file_A;  ///< Distance-function file for anode.
    extern const char* dsF_file_C;  ///< Distance-function file for cathode.

    // -------------------------------------------------------------------------
    // Discretization parameters
    // -------------------------------------------------------------------------

    extern const int    order; ///< Polynomial order for FE basis.
    extern const double dh;    ///< Grid spacing (distance-field scale factor).
    extern const double zeta;  ///< Interface thickness parameter for SBM.
    extern const double ze;    ///< Extrapolation parameter for ψ smoothing.

    // -------------------------------------------------------------------------
    // Numerical tolerances and thresholds
    // -------------------------------------------------------------------------

    extern const double thres; ///< Threshold for phase-field cutoff.
    extern const double eps;   ///< Small epsilon used to avoid division-by-zero.

    // -------------------------------------------------------------------------
    // Time-stepping
    // -------------------------------------------------------------------------

    extern const double dt; ///< Default time step size.
    extern const double tm; ///< Maximum simulation time.

    // -------------------------------------------------------------------------
    // Electrochemical constants
    // -------------------------------------------------------------------------

    extern const double t_minus; ///< Cation transference number.
    extern const double D0;      ///< Base diffusivity.
    extern const double Frd;     ///< Faraday constant scaling factor.
    extern const double Cst1;    ///< Constant used in electrolyte potential transport term.
    extern const double alp;     ///< Charge-transfer coefficient (α).
    extern const double rho_A;   ///< Density/normalization for anode active material.
    extern const double rho_C;   ///< Density/normalization for cathode active material.
    extern const double Cr;      ///< C-rate or reaction scaling constant.
    extern const double gc;      ///< Conductivity or coupling constant (context-dependent).

    // -------------------------------------------------------------------------
    // Voltage / potential limits
    // -------------------------------------------------------------------------

    extern const double Vsr0; ///< Initial cell voltage offset.
    extern const double VCut; ///< Cutoff voltage for stopping criteria.

    // -------------------------------------------------------------------------
    // Initial conditions
    // -------------------------------------------------------------------------

    extern const double init_CnA; ///< Initial anode concentration.
    extern const double init_CnC; ///< Initial cathode concentration.
    extern const double init_CnE; ///< Initial electrolyte concentration.

    extern const double init_BvA; ///< Initial boundary value for anode potential.
    extern const double init_BvC; ///< Initial boundary value for cathode potential.
    extern const double init_BvE; ///< Initial boundary value for electrolyte potential.

    extern const double init_Rxn; ///< Initial reaction rate (global).
    extern const double init_RxA; ///< Initial anode reaction rate.
    extern const double init_RxC; ///< Initial cathode reaction rate.

} // namespace Constants

#endif // CONSTANTS_HPP
