#ifndef CONCENTRATIONS_HPP
#define CONCENTRATIONS_HPP

// Public: Members can be accessed from anywhere. This is the default access modifier. 
// Protected: Members can be accessed within the class and by classes that inherit from that class. 
// Private: Members can only be accessed within the class that defines them.

/**
 * @file Concentrations.hpp
 * @brief Header file for the Concentrations class, which handles the concentration-related operations in battery simulations.
 * 
 * The Concentrations class performs various computations related to concentration fields, including initialization,
 * reaction terms, diffusion calculations, and salt conservation. This class supports both the particle concentration (CnP)
 * and electrolyte concentration (CnE) in the battery simulation.
 */

#include "mfem.hpp"
#include "Mesh_Handler.hpp"
#include <memory>

/**
 * @class Concentrations
 * @brief Class to handle concentration calculations in a battery simulation.
 * 
 * The Concentrations class computes and manages the concentration fields (`CnP` and `CnE`), 
 * as well as related calculations like diffusivity, reaction terms, and salt conservation. 
 * It uses the MFEM library for finite element computations and supports parallelism through MPI.
 */
class Concentrations {
public:

    /**
     * @brief Constructor for the Concentrations class.
     * 
     * Initializes the concentration fields and sets up the mesh and finite element space.
     * 
     * @param pmesh Pointer to the parallel mesh.
     * @param fespace Pointer to the parallel finite element space.
     * @param mh Reference to the MeshHandler object for accessing mesh data.
     */
    Concentrations(mfem::ParMesh* pmesh, mfem::ParFiniteElementSpace* fespace, MeshHandler &mh);
    
    /**
     * @brief Destructor for the Concentrations class.
     */
    virtual ~Concentrations() = default;
    
    MeshHandler &mesh_handler;
    
    mfem::ParGridFunction psi; ///< Potential field for the particle.                    
    mfem::ParGridFunction pse; ///< Potential field for the electrolyte.
    mfem::ParGridFunction *CeT; ///< Temporary grid function for salt concentration calculations.

protected:
    mfem::ParMesh *pmesh;  ///< Pointer to the parallel mesh.
    mfem::ParFiniteElementSpace *fespace; ///< Pointer to the finite element space.
    std::shared_ptr<mfem::HypreParMatrix> Mmat; ///< Shared pointer to the mass matrix.
    std::shared_ptr<mfem::CGSolver> solver; ///< Shared pointer to the conjugate gradient solver.
    mfem::HypreSmoother smoother; ///< Smoother for preconditioning.

    double infx; ///< Reaction current density.
    double CeC = 0.0; ///< Total salt concentration.
    double gCeC = 0.0; ///< Global salt concentration.
    double CeAvg = 0.0; ///< Average salt concentration.
    double Ce0 = 0.001; ///< Initial salt concentration.

    /**
     * @brief Sets the initial concentration values across the domain.
     * 
     * This method initializes the concentration field with a uniform value, typically used
     * to set the initial state of the particle (`CnP`) or electrolyte (`CnE`) concentration in battery simulations.
     * 
     * @param Cn Reference to the ParGridFunction representing the concentration field.
     * @param initial_value The value to set for all elements in the concentration field.
     */
    void SetInitialConcentration(mfem::ParGridFunction &Cn, double initial_value);


    /**
     * @brief Sets up the solver for the mass matrix computation.
     * 
     * This method configures the finite element mass matrix and initializes the solver and preconditioner
     * to solve systems involving particle or electrolyte potentials (`psx`) in the battery simulation.
     * 
     * @param psx Reference to the ParGridFunction representing the potential field.
     * @param Mmat Shared pointer to the mass matrix, which will be constructed and returned.
     * @param m_solver Reference to the conjugate gradient solver to be configured.
     * @param smoother Reference to the preconditioner for the mass matrix.
     */
    void SetUpSolver(mfem::ParGridFunction &psx, std::shared_ptr<mfem::HypreParMatrix> &Mmat, mfem::CGSolver &m_solver, mfem::HypreSmoother &smoother);
    
    /**
     * @brief Applies a Neumann boundary condition by negating a potential field.
     * 
     * This method imposes a Neumann boundary condition by taking the input potential field (`psx`) 
     * and creating its negated version (`PGF`). Neumann boundary conditions represent fluxes 
     * at the boundaries of the simulation domain.
     * 
     * @param psx Reference to the ParGridFunction representing the input potential field.
     * @param PGF Reference to the ParGridFunction where the negated potential field will be stored.
     */
    void ImposeNeumannBC(mfem::ParGridFunction &psx, mfem::ParGridFunction &PGF);


    /**
     * @brief Creates a scaled reaction field based on an input reaction field.
     * 
     * This method generates a new reaction field (`Rx2`) by scaling an existing reaction field (`Rx1`) by a specified value. 
     * In the context of battery simulations, this represents a reaction rate.
     * 
     * @param Rx1 Reference to the ParGridFunction representing the input reaction field.
     * @param Rx2 Reference to the ParGridFunction where the scaled reaction field will be stored.
     * @param value The scaling factor applied to the reaction field.
     */
    void CreateReaction(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, double value);
    
    /**
     * @brief Constructs a force term for the system based on a given field and boundary conditions.
     * 
     * This method creates a linear form (`Fxx`) representing a force term, incorporating contributions 
     * from the domain and optionally from boundary conditions. The force term is built using a grid function 
     * coefficient derived from the input field (`gfc`).
     * 
     * @param gfc Reference to the ParGridFunction representing the input field (e.g., reaction rates or fluxes).
     * @param Fxx Reference to the ParLinearForm where the assembled force term will be stored.
     * @param boundary Array of integers specifying boundary indices for applying boundary conditions.
     * @param m ProductCoefficient used to scale the boundary contributions, representing a Neumann-type boundary condition.
     * @param apply_boundary_conditions Boolean flag indicating whether to include boundary contributions in the force term.
     */
    void ForceTerm(mfem::ParGridFunction &gfc, mfem::ParLinearForm &Fxx, mfem::Array<int> boundary, mfem::ProductCoefficient m, bool apply_boundary_conditions);
   
    /**
     * @brief Calculates the total reaction and the reaction current density.
     * 
     * This method computes the total reaction by integrating the reaction field (`Rx`) over the domain
     * and calculates the reaction current density normalized by the characteristic length of the domain.
     * 
     * @param Rx Reference to the ParGridFunction representing the reaction field (e.g., reaction rates).
     * @param xCrnt Reference to a double that will store the computed total reaction value for the current process.
     */
    void TotalReaction(mfem::ParGridFunction &Rx, double xCrnt);
    
    /**
     * @brief Computes the diffusivity field for the particle or electrolyte.
     * 
     * This method calculates the diffusivity based on the potential field (`psx`) and the concentration field (`Cn`) 
     * for either the particle or the electrolyte. The resulting diffusivity is encapsulated in a 
     * `GridFunctionCoefficient` for use in finite element computations.
     * 
     * @param psx Reference to the ParGridFunction representing the potential field.
     * @param Cn Reference to the ParGridFunction representing the concentration field.
     * @param particle_electrolyte Boolean flag indicating whether to compute diffusivity for the particle (true) or the electrolyte (false).
     * @return A shared pointer to a GridFunctionCoefficient encapsulating the computed diffusivity.
     */
    std::shared_ptr<mfem::GridFunctionCoefficient> Diffusivity(mfem::ParGridFunction &psx, mfem::ParGridFunction &Cn, bool particle_electrolyte );
    
    /**
     * @brief Assembles the stiffness matrix and linear system for diffusion-based problems.
     * 
     * This method constructs the stiffness matrix (`Kmatx`) and the right-hand side vector (`Fxb`) 
     * based on the concentration field (`Cn`) and diffusivity coefficient (`cDx`). The method also 
     * accounts for boundary conditions and the finite element space for diffusion processes in the system.
     * 
     * @param boundary Array of integers representing the boundary degrees of freedom for the system.
     * @param Cn Reference to the ParGridFunction representing the concentration field.
     * @param Fxx Reference to the ParLinearForm representing the right-hand side of the linear system.
     * @param Kmatx Shared pointer to the stiffness matrix that will be assembled.
     * @param X1v Reference to the HypreParVector that holds the solution vector.
     * @param Fxb Reference to the HypreParVector that will hold the boundary contributions to the right-hand side.
     * @param cDx Pointer to the GridFunctionCoefficient representing the diffusivity coefficient.
     */
    void KMatrix(mfem::Array<int> boundary, mfem::ParGridFunction &Cn, mfem::ParLinearForm &Fxx, std::shared_ptr<mfem::HypreParMatrix> &Kmatx, mfem::HypreParVector &X1v, mfem::HypreParVector &Fxb, mfem::GridFunctionCoefficient *cDx);
    
    /**
     * @brief Computes and conserves the salt concentration in the system.
     * 
     * This method computes the total salt concentration (`CeC`) by integrating the product of 
     * the electrolyte concentration field (`Cn`) and the potential field (`psx`). It then adjusts 
     * the concentration field (`Cn`) to ensure salt conservation across the electrolyte.
     * 
     * @param Cn Reference to the ParGridFunction representing the concentration field of the electrolyte.
     * @param psx Reference to the ParGridFunction representing the potential field of the electrolyte.
     */
    void SaltConservation(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);
    
    /**
     * @brief Calculates the degree of lithiation.
     * 
     * This method computes the degree of lithiation by integrating the product of the concentration (`Cn`)
     * and the potential (`psx`) over the domain. The lithiation degree is normalized by the total potential.
     * 
     * @param Cn Reference to the ParGridFunction representing the concentration.
     * @param psx Reference to the ParGridFunction representing the potential.
     */
    void LithiationCalculation(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);

private:
    // MeshHandler &mesh_handler;

    mfem::ParBilinearForm *M = nullptr; ///< Pointer to the mass matrix.
    mfem::ParGridFunction *Ps_gf; ///< Grid function for the potential field.
    mfem::GridFunctionCoefficient *cP; ///< Coefficient for the potential field.

    mfem::Array<int> boundary_dofs; ///< Boundary degrees of freedom.

    mfem::ParGridFunction TmpF; ///< Temporary grid function for intermediate calculations.

    int nE;                                         ///< Number of elements.
    int nC;                                         ///< Number of corners in each element.
    int nV;                                         ///< Number of vertices.

    const mfem::Vector& EVol;                       ///< Element volumes from MeshHandler.
    double gtPsi;                                   ///< Total Psi from MeshHandler.
    double gtPse;                                   ///< Total Pse from MeshHandler.

    mfem::ParGridFunction *Rxx; ///< Grid function for reaction terms.
    mfem::GridFunctionCoefficient *cXx; ///< Coefficient for the reaction term.

    double geCrnt; ///< Global current.

};

#endif // CONCENTRATIONS_HPP
