#ifndef DOMAIN_PARAMETERS_HPP
#define DOMAIN_PARAMETERS_HPP

#include "mfem.hpp"
#include <memory>
#include <string>
#include <vector>

using namespace std;

/**
 * @file Domain_Parameters.hpp
 * @brief Declares the Domain_Parameters class for setting up and managing
 *        domain-specific parameters in the battery simulation.
 *
 * This class is responsible for initializing grid functions (psi, pse, AvP, AvB),
 * interpolating domain parameters depending on the mesh type, computing
 * phase field totals, and calculating the target current used in the simulation.
 */


/**
 * @class Domain_Parameters
 * @brief Class for managing domain-specific parameters in battery simulations.
 *
 * The Domain_Parameters class initializes and stores phase fields (`psi`, `pse`)
 * as well as auxiliary fields (`AvP`, `AvB`). It also computes integral quantities
 * like `gtPsi`, `gtPse`, and the target current `gTrgI` based on the mesh type
 */
class Domain_Parameters {

public:

    /**
     * @brief Construct a new Domain_Parameters object.
     * 
     * @param geo Reference to the geometry initialization object, used to access
     *        mesh, finite element space, and distance function data.
     */
    Domain_Parameters(Initialize_Geometry &geo);

    Initialize_Geometry &geometry;
    virtual ~Domain_Parameters();

    /**
     * @brief Setup all domain parameters based on the mesh type.
     * 
     * Initializes grid functions, interpolates parameters (psi, pse, AvP, AvB),
     * calculates totals, and computes the target current.
     * 
     * @param mesh_type Character flag for mesh geometry type:
     *        - "r" rectangle
     *        - "c" circle
     *        - "d" disk
     *        - "v" voxel
     */
    void SetupDomainParameters(const char* mesh_type);


    std::unique_ptr<mfem::ParGridFunction> psi; ///< Solid phase potential
    std::unique_ptr<mfem::ParGridFunction> pse; ///< Electrolyte phase potential
    std::unique_ptr<mfem::ParGridFunction> psA; ///< Anode phase potential
    std::unique_ptr<mfem::ParGridFunction> psC; ///< Cathode phase potential

    std::unique_ptr<mfem::ParGridFunction> AvP; ///< Particle surface area
    std::unique_ptr<mfem::ParGridFunction> AvA; ///< Anode surface area
    std::unique_ptr<mfem::ParGridFunction> AvC; ///< Cathode surface area
    std::unique_ptr<mfem::ParGridFunction> AvB; ///< Boundary surface area

    double gtPsi; ///< Global total for Psi
    double gtPse; ///< Global total for Pse
    double gTrgI; ///< Global target current
    mfem::Vector EVol; ///< Element volumes



private:

    void InitializeGridFunctions();

    /**
     * @brief Interpolate the distance function into domain parameters.
     * 
     * Projects psi and AvP depending on mesh type ("r", "c", "d", "v") 
     * and clamps them to avoid unphysical values.
     * 
     * @param mesh_type Mesh type specifier.
     */
    void InterpolateDomainParameters(const char* mesh_type);

    /**
     * @brief Generic routine to compute total integral of a field.
     * 
     * @param grid_function The field to integrate.
     * @param element_volumes Element volumes for weighting.
     * @param local_total Accumulated local total (output).
     * @param global_total Accumulated global total after MPI reduction (output).
     */
    void CalculateTotals(const mfem::ParGridFunction& grid_function, const mfem::Vector& element_volumes, double& local_total, double& global_total);
    
    
    /**
     * @brief Build element volumes then integrate a phase field.
     *
     * Fills @c EVol from the mesh, then calls @ref CalculateTotals to compute
     * the local and global totals for the provided phase field.
     *
     * @param grid_function Phase indicator field (e.g., ψ or ψ_e).
     * @param total         [out] Process-local total.
     * @param global_total  [out] MPI-reduced total across all ranks.
     */
    void CalculateTotalPhaseField(const mfem::ParGridFunction& grid_function, double& total, double& global_total);
    
    /**
     * @brief Compute global ψ and ψ_e integrals and derive the target current.
     *
     * Calls @ref CalculateTotalPhaseField for both ψ and ψ_e, then updates
     * the global target current via @ref CalculateTargetCurrent.
     */
    void CalculatePhasePotentialsAndTargetCurrent();
    
    /**
     * @brief Derive the (local) target current from integrated ψ.
     *
     * Uses model constants (ρ, Cr, etc.) to compute a current target from the
     * integrated particle phase volume and then participates in a global reduction.
     *
     * @param total_psi Local integral of ψ before MPI reduction.
     */
    void CalculateTargetCurrent(double total_psi);
    
    /**
     * @brief Print a summary of domain totals (rank 0 only).
     *
     * Logs gtPsi, gtPse, and gTrgI to stdout for quick inspection.
     */
    void PrintInfo();
    
    int nV; ///< Number of vertices in the mesh
    int nE; ///< Number of elements in the mesh
    int nC; ///< Number of corners per element
    
    mfem::ParGridFunction* dsF;   ///< Pointer to distance function grid
    mfem::ParGridFunction* dsF_A; ///< Pointer to anode distance function grid
    mfem::ParGridFunction* dsF_C; ///< Pointer to cathode distance function grid
    mfem::ParMesh* pmesh;         ///< Pointer to parallel mesh
    std::shared_ptr<mfem::ParFiniteElementSpace> fespace; ///< Shared pointer to finite element space


    double tPsi; ///< Target Psi value
    double tPse; ///< Target Pse value
    double trgI; ///< Target current


};




#endif //DOMAIN_PARAMETERS_HPP