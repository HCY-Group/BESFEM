// PotentialBase.hpp
#pragma once
#include "mfem.hpp"
#include "Initialize_Geometry.hpp"
#include "Domain_Parameters.hpp"
#include "BoundaryConditions.hpp"
#include "Utils.hpp"
#include "FEMOperators.hpp"
#include <memory>

class PotentialBase {
public:
    
    PotentialBase(Initialize_Geometry &geo, Domain_Parameters &para);

    virtual ~PotentialBase() = default;

    virtual void SetupField(mfem::ParGridFunction &ph, double initial_value, mfem::ParGridFunction &psx) = 0;
    virtual void AssembleSystem(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &potential) = 0;
    virtual void UpdatePotential(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror) = 0;

    // double Xfr = 0.0;
    // double CeC = 0.0;
    // double gCeC = 0.0;                                                ///< Global total salt
    // mfem::Array<double> VtxVal;                                        ///< Vertex values buffer.


    int    nE = 0;                            ///< Number of elements.
    int    nC = 0;                            ///< Nodal count per element (corners).
    int    nV = 0;                            ///< Number of vertices (global).
    // mfem::Vector EAvg;                        ///< Per-element averages (work buffer).
    const mfem::Vector &EVol;                 ///< Element volumes (from domain parameters).
    // double CeAvg = 0.0;                                               ///< Global average salt.
    // double Ce0   = 0.001;                                             ///< Initial salt reference.

    // double gtPsi = 0.0;                                                ///< Global ψ_P normalization.
    // double gtPse = 0.0;                                                ///< Global ψ_E normalization.
    // mfem::Array<int> boundary_dofs; ///< Boundary degrees of freedom
    // mfem::HypreParVector X1v; ///< Scratch vector

    Initialize_Geometry &geometry;            ///< Geometry/mesh context
    Domain_Parameters   &domain_parameters;   ///< Domain/physics parameters

    mfem::ParMesh *pmesh = nullptr;                                  ///< Parallel mesh
    // mfem::Mesh    *gmesh = nullptr;                                   ///< Serial/global mesh

    std::shared_ptr<mfem::ParFiniteElementSpace> fespace;             ///< Parallel FE space
    mfem::ParGridFunction TmpF; ///< Temporary field for error calculations

    mfem::ParGridFunction px0; ///< Pointer to the potential grid function
    mfem::HypreParVector X0; ///< Solution vector



};





















// #ifndef POTENTIALS_HPP
// #define POTENTIALS_HPP

// // Public: Members can be accessed from anywhere. This is the default access modifier.
// // Protected: Members can be accessed within the class and by classes that inherit from that class. 
// // Private: Members can only be accessed within the class that defines them.

// /**
//  * @file Potentials.hpp
//  * @brief Base class for assembling/solving potential fields in the battery simulation.
//  *
//  * Provides common helpers to initialize potentials, assemble forcing terms,
//  * compute global errors, and access shared mesh/space context.
//  */

// #include "mfem.hpp"
// #include "Initialize_Geometry.hpp"
// // #include "SolverSteps.hpp"
// #include "Domain_Parameters.hpp"
// #include "BoundaryConditions.hpp"
// #include "FEMOperators.hpp"

// #include <memory>

// /**
//  * @class Potentials
//  * @brief Base class for potential-field computations (electrode/electrolyte).
//  *
//  * Inherits low-level assembly/solver utilities from @ref SolverSteps and
//  * exposes helpers shared by concrete potential solvers (e.g., PotE, PotP).
//  */
// class Potentials {
// public:

//     /**
//      * @brief Construct the potentials helper.
//      * @param geo  Geometry/space container (mesh, FE spaces, counts).
//      * @param para Domain/physics parameters (constants, totals, etc.).
//      */
//     Potentials(Initialize_Geometry &geo, Domain_Parameters &para);

//     virtual ~Potentials() = default; ///< Virtual destructor

//     Initialize_Geometry &geometry;    ///< Geometry/mesh context
//     Domain_Parameters &domain_parameters;  ///< Domain/physics parameters

//     /**
//      * @brief Initialize a potential field to a uniform value.
//      * @param ph            Potential grid function (in/out).
//      * @param initial_value Value assigned to all DoFs.
//      */
//     void SetInitialPotentials(mfem::ParGridFunction &ph, double initial_value);


//     /**
//      * @brief Assemble a force (RHS) vector from a scaled reaction field.
//      *
//      * Forms domain/boundary contributions using a coefficient built from
//      * the reaction field and writes into the provided RHS forms.
//      *
//      * @param Rx1       Input reaction field.
//      * @param Rx2       Output scaled reaction field.
//      * @param value     Scale factor applied to @p Rx1.
//      * @param coef      Coefficient wrapping the (possibly updated) field.
//      * @param rhs_form  Owned linear form to (re)assemble.
//      * @param rhs_form2 Target linear form to receive assembled values.
//      */
//     void AssembleForceVector(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, double value, mfem::GridFunctionCoefficient &coef, std::unique_ptr<mfem::ParLinearForm> &rhs_form, mfem::ParLinearForm &rhs_form2);
    

    
//     /**
//      * @brief Compute a global error metric for a potential field.
//      *
//      * Typically integrates a weighted norm of (reference - current) using
//      * the phase mask and element volumes; reduces across MPI.
//      *
//      * @param px0         Reference/previous potential (input).
//      * @param potential   Current potential (input).
//      * @param psx         Phase mask used for weighting (input).
//      * @param globalerror [out] MPI-reduced error value.
//      * @param gtPsx       Global normalization factor for weighting.
//      */
//     void ComputeGlobalError(mfem::ParGridFunction &px0, mfem::ParGridFunction &potential, mfem::ParGridFunction &psx, double &globalerror, double gtPsx);

//     double Vcell; ///< Cell voltage

//     int nE; ///< Number of elements in the mesh
//     int nC; ///< Number of nodes per element (corners)
//     int nV; ///< Number of vertices in the mesh

// protected:
    
//     mfem::ParMesh *pmesh; ///< Pointer to the parallel mesh
//     std::shared_ptr<mfem::ParFiniteElementSpace> fespace; ///< Pointer to the finite element space


// private:

//     std::unique_ptr<mfem::ParGridFunction> Rxx; ///< Intermediate grid function for reactions

//     std::unique_ptr<mfem::ParLinearForm> Bx2; ///< Scratch RHS form.
//     std::unique_ptr<mfem::GridFunctionCoefficient> cXx; ///< Coefficient derived from grid functions

//     mfem::ParGridFunction px0; ///< Pointer to the potential grid function
//     mfem::HypreParVector X0; ///< Solution vector


//     mfem::ParGridFunction TmpF; ///< Temporary field for error calculations

//     const mfem::Vector& EVol; ///< Reference to element volumes from the MeshHandler

// };

// #endif // POTENTIALS_HPP
