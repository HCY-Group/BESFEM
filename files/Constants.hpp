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
    
    /**
     * @brief Path to the mesh file used in the simulation.
     */
    extern const char* mesh_file;

    /**
     * @brief Path to the file containing data for the surface flux (dsF).
     */
    extern const char* dsF_file;

    /**
     * @brief Order of the finite element basis functions.
     */
    extern const int order;

    /**
     * @brief Mesh element size.
     */
    extern const double dh;

    /**
     * @brief Interfacial thickness.
     */
    extern const double zeta;

    /**
     * @brief Threshold value used in various numerical checks.
     */
    extern const double thres;

    /**
     * @brief Small epsilon value for numerical tolerance.
     */
    extern const double eps;
    
    /**
     * @brief Time step size.
     */
    extern const double dt;

    /**
     * @brief Initial simulation time.
     */
    extern const double tm;

    /**
     * @brief Transference number (ratio of cation conductivity to total conductivity).
     */
    extern const double t_minus;

    /**
     * @brief Base diffusivity.
     */
    extern const double D0;

    /**
     * @brief Faraday constant in C/mol (electric charge per mole of electrons).
     */
    extern const double Frd;

    /**
     * @brief Constant for exponential calculations in the Nernst equation.
     */
    extern const double Cst1;

    /**
     * @brief Symmetry factor for charge transfer kinetics.
     */
    extern const double alp;

    /**
     * @brief Lithium site density.
     */
    extern const double rho;

    /**
     * @brief C-rate (charging/discharging rate normalized to battery capacity).
     */
    extern const double Cr;

    /**
     * @brief Voltage scanning rate.
     */
    extern const double Vsr;

    /**
     * @brief Cut-off voltage for simulation termination.
     */
    extern const double VCut;
}

#endif // CONSTANTS_HPP
