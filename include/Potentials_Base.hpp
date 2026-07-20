/**
 * @file Potentials_Base.hpp
 * @brief Defines the abstract base class for electric-potential solvers in BESFEM.
 */

#pragma once

#include "mfem.hpp"
#include "Initialize_Geometry.hpp"
#include "Domain_Parameters.hpp"
#include "BoundaryConditions.hpp"
#include "Utils.hpp"
#include "FEMOperators.hpp"
#include <memory>
#include <vector>

/**
 * @class PotentialBase
 * @brief Abstract base class for BESFEM electric-potential solvers.
 *
 * PotentialBase defines the shared interface and common FEM data used by
 * electrode and electrolyte potential solvers. Derived classes initialize
 * potential fields, assemble material-dependent systems, solve the finite
 * element equations, and update boundary voltages when needed.
 */
class PotentialBase {
public:

    /**
     * @brief Construct the base potential solver.
     *
     * Stores references to the geometry and domain-parameter objects and
     * initializes shared mesh, finite-element-space, and element-volume data.
     *
     * @param geo Reference to the initialized geometry object.
     * @param para Reference to the domain-parameter object.
     */
    PotentialBase(Initialize_Geometry &geo, Domain_Parameters &para);

    /// Virtual destructor for safe cleanup through base-class pointers.
    virtual ~PotentialBase() = default;

    /**
     * @brief Initialize the potential field and solver-specific operators.
     *
     * Derived classes implement boundary-condition setup, initial potential
     * assignment, phase-field masking, and any matrix/preconditioner
     * initialization needed before timestepping.
     *
     * @param ph Potential field to initialize.
     * @param initial_value Scalar initial potential.
     * @param psx Phase-field mask for the active region.
     */
    virtual void SetupField(mfem::ParGridFunction &ph,
                            double initial_value,
                            mfem::ParGridFunction &psx) = 0;

    /**
     * @brief Assemble the potential system for a single concentration field.
     *
     * Derived classes build the finite-element matrix and right-hand side using
     * the current concentration, phase-field mask, and potential field.
     *
     * @param Cn Concentration field used to evaluate material properties.
     * @param psx Phase-field mask for the active region.
     * @param potential Potential field associated with the solve.
     */
    virtual void AssembleSystem(mfem::ParGridFunction &Cn,
                                mfem::ParGridFunction &psx,
                                mfem::ParGridFunction &potential) = 0;

    /**
     * @brief Assemble the potential system for multiple particle/material groups.
     *
     * The default implementation may be overridden by derived classes that
     * support multi-particle or multi-material conductivity assembly.
     *
     * @param Cn_groups Concentration fields for each particle/material group.
     * @param psi_groups Phase-field masks for each particle/material group.
     * @param materials Material type associated with each group.
     * @param potential Potential field associated with the solve.
     */
    virtual void AssembleSystem(
        const std::vector<mfem::ParGridFunction*> &Cn_groups,
        const std::vector<mfem::ParGridFunction*> &psi_groups,
        const std::vector<sim::MaterialType> &materials,
        mfem::ParGridFunction &potential);

    /**
     * @brief Solve for the updated potential field.
     *
     * Derived classes solve the assembled potential equation, update the
     * potential field, and compute a global error metric.
     *
     * @param Rx Reaction/source field.
     * @param phx Potential field to update.
     * @param psx Phase-field mask for the active region.
     * @param gerror Output global error metric.
     */
    virtual void UpdatePotential(mfem::ParGridFunction &Rx,
                                 mfem::ParGridFunction &phx,
                                 mfem::ParGridFunction &psx,
                                 double &gerror) = 0;

    /**
     * @brief Return the current applied boundary voltage.
     *
     * @return Boundary voltage used by the potential solver.
     */
    virtual double GetBoundaryVoltage() const = 0;

    /**
     * @brief Increment the applied boundary voltage.
     *
     * @param dV Voltage increment to add.
     */
    virtual void AddBoundaryVoltage(double dV) = 0;

    virtual mfem::ParGridFunction GetConductivity() = 0;

    int nE = 0; ///< Number of mesh elements.
    int nC = 0; ///< Number of nodes per element.
    int nV = 0; ///< Number of vertices.

    const mfem::Vector &EVol; ///< Element volumes from Domain_Parameters.

    Initialize_Geometry &geometry;          ///< Geometry and mesh data.
    Domain_Parameters   &domain_parameters; ///< Domain fields and global quantities.

    mfem::ParMesh *pmesh = nullptr; ///< Parallel mesh pointer.

    std::shared_ptr<mfem::ParFiniteElementSpace> fespace; ///< Parallel finite element space.

    mfem::ParGridFunction TmpF; ///< Temporary grid-function workspace.

    mfem::ParGridFunction px0; ///< Previous or working potential field.
    mfem::HypreParVector  X0;  ///< True-DOF vector for the potential solve.
};
