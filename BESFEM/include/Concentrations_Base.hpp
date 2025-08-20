#ifndef CONCENTRATIONS_HPP
#define CONCENTRATIONS_HPP

// Public: Members can be accessed from anywhere. This is the default access modifier. 
// Protected: Members can be accessed within the class and by classes that inherit from that class. 
// Private: Members can only be accessed within the class that defines them.

/**
 * @file Concentrations.hpp
 * @brief Header file for the Concentrations class, which handles the concentration-related operations in battery simulations.
 */

#include "mfem.hpp"
#include "Initialize_Geometry.hpp"
#include "SolverSteps.hpp"
#include "Domain_Parameters.hpp"
#include <memory>

/**
 * @defgroup concentrations Concentration Modules
 * @brief Classes that advance concentration fields.
 * @{
 */

/**
 * @class Concentrations
 * @brief Base class for concentration calculations (shared by CnCH/CnE/etc.).
 * @ingroup concentrations
 *
 * @details
 * - Holds shared geometry/space, buffers, and coefficients.
 * - Implements helpers for initialization, reaction assembly, diffusivity,
 *   salt conservation, and lithiation computation.
 * - Inherits low-level solver utilities from @ref SolverSteps.
 */
class Concentrations : public SolverSteps {
public:

   /**
     * @brief Construct the base concentration helper.
     * @param geo  Geometry/space container (mesh, FE space, counts).
     * @param para Domain/physics parameters (volumes, normalizations, etc.).
     */
    Concentrations(Initialize_Geometry &geo, Domain_Parameters &para);

    double GetLithiation() const { return Xfr; } ///< Return current degree of lithiation.

    
    // ---------- Geometry / spaces ----------
    Initialize_Geometry &geometry;            ///< Geometry/mesh context
    Domain_Parameters   &domain_parameters;   ///< Domain/physics parameters

    
    // ---------- Fields / masks ----------
    mfem::ParGridFunction psi;                ///< Particle phase mask/weight ψ_P.
    mfem::ParGridFunction pse;                ///< Electrolyte phase mask/weight ψ_E.
    mfem::ParGridFunction *CeT = nullptr;     ///< (Optional) temp field for salt conservation.

    mfem::Array<int> boundary_dofs; ///< Boundary degrees of freedom

    int    nE = 0;                            ///< Number of elements.
    int    nC = 0;                            ///< Nodal count per element (corners).
    int    nV = 0;                            ///< Number of vertices (global).
    double Xfr = 0.0;                         ///< Degree of lithiation.
    double infx = 0.0;                        ///< Reaction current density (flux scale).
    mfem::Vector EAvg;                        ///< Per-element averages (work buffer).
    const mfem::Vector &EVol;                 ///< Element volumes (from domain parameters).
    double geCrnt = 0.0;                      ///< Global reaction current (reduced).


    /**
     * @brief Enforce salt conservation across the electrolyte.
     * @param Cn  Electrolyte concentration grid function (in/out).
     * @param psx Phase field ψ_E used for masking.
     */
    void SaltConservation(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);


protected:
   // ---------- Mesh / FE / solver  ----------
    mfem::ParMesh *pmesh = nullptr;                                  ///< Parallel mesh
    mfem::Mesh    *gmesh = nullptr;                                   ///< Serial/global mesh
    std::shared_ptr<mfem::ParFiniteElementSpace> fespace;             ///< Parallel FE space

    std::shared_ptr<mfem::HypreParMatrix> Mmat;                       ///< Mass matrix 
    std::shared_ptr<mfem::CGSolver>       solver;                     ///< CG solver 
    mfem::HypreSmoother                   smoother;                   ///< Preconditioner smoother.

    double CeC  = 0.0;                                                ///< Local total salt.
    double gCeC = 0.0;                                                ///< Global total salt
    double CeAvg = 0.0;                                               ///< Global average salt.
    double Ce0   = 0.001;                                             ///< Initial salt reference.
    double L_w   = 0.0;                                               ///< Characteristic boundary size for flux.

    /**
     * @brief Initialize concentration field uniformly.
     * @param Cn            Concentration field (in/out).
     * @param initial_value Value to assign to all DoFs.
     */
    void SetInitialConcentration(mfem::ParGridFunction &Cn, double initial_value);

    
    /**
     * @brief Impose Neumann-type condition by negating a field.
     * @param psx Input phase/potential field.
     * @param PGF Output negated field.
     */
    void ImposeNeumannBC(mfem::ParGridFunction &psx, mfem::ParGridFunction &PGF);


    /**
     * @brief Create a scaled reaction field Rx2 = value * Rx1.
     * @param Rx1  Input reaction field.
     * @param Rx2  Output scaled reaction field.
     * @param value Scale factor.
     */
    void CreateReaction(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, double value);
    
    /**
     * @brief Compute total reaction and flux scale (infx).
     * @param Rx    Reaction field (in).
     * @param xCrnt Local process reaction total (out).
     */
    void TotalReaction(mfem::ParGridFunction &Rx, double &xCrnt);
    
    /**
     * @brief Build a diffusivity coefficient from Cn and ψ.
     * @param psx                 Phase mask ψ.
     * @param Cn                  Concentration field.
     * @param particle_electrolyte If true, use particle model; else electrolyte.
     * @return Shared coefficient wrapping the computed diffusivity field.
     */
    std::shared_ptr<mfem::GridFunctionCoefficient> Diffusivity(mfem::ParGridFunction &psx, mfem::ParGridFunction &Cn, bool particle_electrolyte );

    
    /**
     * @brief Compute the lithiation degree Xfr from Cn and ψ.
     * @param Cn  Concentration field.
     * @param psx Phase field mask.
     */
    void LithiationCalculation(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);

private:
    // ---------- Internal assembly helpers / storage ----------
    mfem::ParBilinearForm *M = nullptr;                               ///< Mass bilinear form (owned elsewhere).
    mfem::ParGridFunction  *Ps_gf = nullptr;                           ///< Potential/phase grid function handle.
    mfem::GridFunctionCoefficient *cP = nullptr;                       ///< Coefficient wrapping Ps_gf.

    std::unique_ptr<mfem::ParGridFunction> TmpF;                       ///< Temporary field for products.

    double gtPsi = 0.0;                                                ///< Global ψ_P normalization.
    double gtPse = 0.0;                                                ///< Global ψ_E normalization.

    mfem::ParGridFunction *Rxx = nullptr;                              ///< Reaction field handle.
    mfem::GridFunctionCoefficient *cXx = nullptr;                      ///< Coefficient wrapping reaction.
    double inv_nC = 0.0;                                               ///< 1.0 / nC (cached).

    mfem::Array<double> VtxVal;                                        ///< Vertex values buffer.
    std::unique_ptr<mfem::ParLinearForm> Bx2;                          ///< Linear form (aux).
    std::shared_ptr<mfem::GridFunctionCoefficient> last_cDx;           ///< Last computed diffusivity coeff.
    std::shared_ptr<mfem::ParGridFunction> Dx;                          ///< Persistent diffusivity field.

};

/** @} */ // end group concentrations


#endif // CONCENTRATIONS_HPP
