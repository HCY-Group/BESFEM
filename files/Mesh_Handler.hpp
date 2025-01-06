#ifndef MESH_HANDLER_HPP
#define MESH_HANDLER_HPP

/**
 * @file Mesh_Handler.hpp
 * @brief Header file for the MeshHandler class, which manages mesh operations, finite element spaces, and related computations.
 * 
 * This file defines the MeshHandler class and its methods for handling mesh initialization, grid functions,
 * and domain parameter interpolation using the MFEM library.
 */


#include "Constants.hpp"
#include "mfem.hpp"
#include <memory>

using namespace mfem;

extern double gTrgI;

/**
 * @class MeshHandler
 * @brief Class for managing mesh operations and domain-specific computations.
 * 
 * The MeshHandler class initializes and processes finite element meshes,
 * manages grid functions, and computes domain parameters for parallel simulations.
 */
class MeshHandler {

public:
    /**
     * @brief Default constructor
     * 
     * Initializes simulation parameters and file paths using the Constants module
     */
    MeshHandler();

    /**
     * @brief Loads and initializes the mesh
     * 
     * Combines mesh initialization, mesh information printing, and MPI synchronization
     */
    void LoadMesh();

    /**
     * @brief Sets up boundary conditions for the simulation
     * 
     * @param[in] pmesh Pointer to the parallel mesh object
     * @param[in] fespace Pointer to the parallel finite element space object
     */
    void SetupBoundaryConditions(mfem::ParMesh *pmesh, mfem::ParFiniteElementSpace *fespace);

    std::unique_ptr<mfem::ParGridFunction> psi; ///< Solid phase potential
    std::unique_ptr<mfem::ParGridFunction> pse; ///< Electrolyte phase potential
    std::unique_ptr<mfem::ParGridFunction> AvP; ///< Particle surface area
    std::unique_ptr<mfem::ParGridFunction> AvB; ///< Boundary surface area

    int nV; ///< Number of vertices
    int nE; ///< Number of elements
    int nC; ///< Number of corners per element

    double L_w; ///< West boundary size

    mfem::Array<int> nbc_w_bdr; ///< West Neumann Boundary Conditions
    mfem::Array<int> ess_tdof_list_w; ///< Total DOF West
    mfem::Array<int> ess_tdof_list_e; ///< Total DOF East

    mfem::Array<int> dbc_w_bdr; ///< West Dirichlet Boundary Conditions
    mfem::Array<int> dbc_e_bdr; ///< East Dirichlet Boundary Conditions

    std::unique_ptr<mfem::ParMesh> pmesh; ///< Parallel mesh

    mfem::Vector EVol; ///< Element volumes

    double gtPsi; ///< Global total for Psi
    double gtPse; ///< Global total for Pse


private:
    // Member functions

    /**
     * @brief Initializes the mesh and associated data structures
     * 
     * Sets up global and local meshes, calculates element volumes, initializes grid functions,
     * and interpolates domain-specific parameters
     */
    void InitializeMesh();

    /**
     * @brief Calculates the volume of each element in the mesh
     * 
     * @param[in] pmesh Pointer to the parallel mesh
     * @param[out] EVol Vector to store calculated element volumes
     */
    void CalculateElementVolume(const std::unique_ptr<ParMesh>& pmesh, Vector& EVol);

    /**
     * @brief Reads the global distance function from a file and stores it in a grid function
     * 
     * @param[in] fespace Pointer to the finite element space used for the global grid function
     */
    void ReadGlobalDistanceFunction(const std::unique_ptr<FiniteElementSpace>& fespace);

    /**
     * @brief Initializes grid functions for the domain simulation
     * 
     * @param[in] fespace Pointer to the parallel finite element space
     */
    void InitializeGridFunctions(const std::shared_ptr<mfem::ParFiniteElementSpace>& fespace);

    /**
     * @brief Interpolates domain parameters using the distance function
     * 
     * Updates Psi, Pse, and their derivatives based on domain-specific parameters
     * @param[in] fespace Pointer to the parallel finite element space
     */
    void InterpolateDomainParameters(const std::shared_ptr<ParFiniteElementSpace>& fespace);

    /**
     * @brief Calculates totals for a given grid function
     * 
     * @param[in] grid_function Reference to the grid function for which totals are calculated
     * @param[in] element_volumes Vector containing element volumes
     * @param[out] local_total Total value calculated locally on the process
     * @param[out] global_total Total value aggregated across all processes
     */
    void CalculateTotals(const mfem::ParGridFunction& grid_function, const mfem::Vector& element_volumes, double& local_total, double& global_total);
    
    /**
     * @brief Calculates the total value of a phase field (e.g., Psi or Pse) across the domain
     * 
     * This method computes the local and global totals for a given grid function
     * by integrating its values over the mesh elements
     * 
     * @param[in] grid_function Reference to the grid function representing the phase field
     * @param[out] total Local total value computed on the current process
     * @param[out] global_total Global total value aggregated across all processes using MPI
     */
    void CalculateTotalPhaseField(const mfem::ParGridFunction& grid_function, double& total, double& global_total);
    
    /**
     * @brief Calculates the total Psi and Pse fields and computes the target current
     */
    void CalculatePhasePotentialsAndTargetCurrent();
    
    /**
     * @brief Computes the target current based on the total Psi value
     * 
     * @param[in] total_psi The total Psi value computed across the domain
     */
    void CalculateTargetCurrent(double total_psi);
    
    /**
     * @brief Prints mesh-related information, including size and calculated values
     */
    void PrintMeshInfo();

    // Member variables
    const char* mesh_file; ///< Path to the mesh file
    const char* dsF_file; ///< Path to the distance function file
    int order; ///< Polynomial order of the finite element basis functions
    
    double dh; ///< Mesh element size from Constants
    double zeta; ///< Interfacial thickness from Constants
    double eps; ///< Epsilon value from Constants
    double rho; ///< Lithium site density from Constants
    double Cr; ///< C-rate from Constants
    
    double Onm; ///< Number of grid function entries
    double localTotal; ///< Total computed value
    double val; ///< Intermediate value
    double globalTotal; ///< Global sum
    double tPsi; ///< Target Psi value
    double tPse; ///< Target Pse value
    double trgI; ///< Target current

    // Mesh and FE space
    std::unique_ptr<mfem::Mesh> gmesh; ///< initial global mesh
    std::unique_ptr<mfem::FiniteElementSpace> gFespace; ///< global finite element space
    std::shared_ptr<mfem::ParFiniteElementSpace> fespace; ///< finite element space
    std::unique_ptr<mfem::ParMesh> pmesh0; ///< initial local pmesh
    std::unique_ptr<mfem::ParGridFunction> dsF; ///< local distance function
    std::unique_ptr<mfem::GridFunction> gDsF; ///< Global distance function grid

};

#endif // MESH_HANDLER_HPP
