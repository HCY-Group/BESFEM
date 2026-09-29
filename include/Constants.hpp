#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

/**
 * @file Constants.hpp
 * @brief Declares global constants used throughout BESFEM simulations.
 *
 * Includes polynomial order, numerical thresholds, electrochemical constants,
 * and initial reaction values. Runtime settings are stored in SimulationConfig.
 */

/**
 * @namespace Constants
 * @brief Compiled numerical and electrochemical constants.
 *
 * This namespace provides a central location for:
 * - finite element order,
 * - physical parameters,
 * - model constants,
 * - default initialization values for reactions.
 *
 * SimulationConfig provides separate runtime settings and uses order as a default.
 */
namespace Constants {

    // -------------------------------------------------------------------------
    // Discretization parameters
    // -------------------------------------------------------------------------

    extern const int    order; ///< Polynomial order for FE basis.
    extern const double zeta;  ///< Interface thickness parameter for SBM.

    // -------------------------------------------------------------------------
    // Numerical tolerances and thresholds
    // -------------------------------------------------------------------------

    extern const double thres; ///< Threshold for phase-field cutoff.
    extern const double eps;   ///< Small epsilon used to avoid division-by-zero.

    // -------------------------------------------------------------------------
    // Electrochemical constants
    // -------------------------------------------------------------------------

    extern const double t_minus; ///< Cation transference number.
    extern const double D0;      ///< Base diffusivity.
    extern const double Frd;     ///< Faraday constant scaling factor.
    extern const double Cst1;    ///< Constant used in electrolyte potential transport term.
    extern const double alp;     ///< Charge-transfer coefficient (α).

    // -------------------------------------------------------------------------
    // Initial conditions
    // -------------------------------------------------------------------------

    extern const double init_Rxn; ///< Initial reaction rate (global).
    extern const double init_RxA; ///< Initial anode reaction rate.
    extern const double init_RxC; ///< Initial cathode reaction rate.

    extern const double RT;   ///< RT constant at 300K.
    extern const double Perm; ///< Permittivity constant.

} // namespace Constants

#endif // CONSTANTS_HPP
