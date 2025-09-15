#ifndef POTA_HPP
#define POTA_HPP

/**
 * @file PotA.hpp
 * @brief Derived class for anode potential solver in the battery simulation.
 *
 * Implements initialization, conductivity assembly, and time-stepping for
 * the particle (solid-phase) potential field. Extends the base Potentials
 * class to specialize for electrode domain physics.
 */

#include "Potentials_Base.hpp"

/**
 * @brief Boundary value for the particle potential
 */
extern double BvA;

/**
 * @defgroup potentials Potential Modules
 * @brief Classes that assemble and advance potential fields.
 * @{
 */

/**
 * @class PotA
 * @brief Derived class implementing the particle potential model for battery simulations
 * @ingroup potentials
 *
 * This class extends the `Potentials` base class to handle computations
 * related to particle potentials, including initialization, time-stepping,
 * and global error calculations
 */
class PotA : public Potentials {



public:

    /**
     * @brief Constructor for the PotP class
     * @param pm Pointer to the parallel mesh
     * @param fe Pointer to the finite element space
     * @param mh Reference to the mesh handler
     */
    PotA(Initialize_Geometry &geo, Domain_Parameters &para);


    Initialize_Geometry &geometry; ///< Geometry/mesh handler
    Domain_Parameters &domain_parameters;  ///< Domain parameters

    double BvA; ///< Boundary value for particle potential


    /**
     * @brief Initializes the particle potential field and solver
     * @param ph Particle potential grid function
     * @param initial_value Initial potential value
     */
    void Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx);

    /**
     * @brief Performs time-stepping for the particle potential
     * @param Cn Particle concentration field
     * @param psx Psi potential field
     * @param phx Particle potential field
     */
    void TimeStep(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &potential);
    
    
    /**
     * @brief Advance particle potential by one time step.
     *
     * @param Cn        Particle concentration field.
     * @param psx       Psi phase-field mask.
     * @param potential Particle potential field (solution, in/out).
     */
    void Advance(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror);



private:

    void ParticleConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx); ///< Computes particle conductivity

    std::shared_ptr<mfem::ParFiniteElementSpace> fespace; ///< Pointer to the finite element space

    mfem::CGSolver cgPP_solver; ///< Conjugate gradient solver for the particle potential system

    mfem::ParLinearForm B1t; ///< Linear form for the right-hand side
    mfem::HypreParVector X1v; ///< Solution vector
    mfem::HypreParVector B1v; ///< Right-hand-side vector
    mfem::HypreParVector Fpb; ///< Right-hand-side vector
    mfem::HypreParVector Xs0; ///< Solution vector for particle potential
    std::unique_ptr<mfem::HypreBoomerAMG> Mpp; ///< BoomerAMG preconditioner (built once with KmP)



    double gtPsi; ///< Total Psi from MeshHandler

    mfem::Array<int> ess_tdof_list_e; ///< List of essential true degrees of freedom for Dirichlet boundary conditions
    mfem::Array<int> dbc_e_bdr; ///< Array marking Dirichlet boundary attributes

    mfem::ParGridFunction kap; ///< Conductivity field for particle potential
    mfem::GridFunctionCoefficient cKp; ///< Coefficient for the conductivity field
    mfem::GridFunctionCoefficient cRp; ///< Coefficient for the reaction field
    mfem::ParGridFunction RpP; ///< Reaction field for the particle potential
    mfem::ParGridFunction pP0; ///< Particle potential grid function

    std::unique_ptr<mfem::ParBilinearForm> Kp2; ///< Bilinear form for particle potential conductivity
    std::unique_ptr<mfem::ParLinearForm> Bp2; ///< Linear form for the reaction term
    mfem::HypreParMatrix KmP; ///< Stiffness matrix for particle potential conductivity

    mfem::ParLinearForm Fpt; ///< Linear form for the force term in particle potential calculations

};

/** @} */ // end group potentials


#endif // POTA_HPP