// CnCH.hpp
#ifndef CNCH_HPP
#define CNCH_HPP

#include "mfem.hpp"
#include "Concentrations_Base.hpp"
// #include "../code/Initialize_Geometry.hpp"
// #include "..Domain_Parameters.hpp"
// #include "../code/SolverSteps.hpp"
#include <memory>
#include <vector>

/**
 * @class CnCH
 * @brief Implements Cahn-Hilliard concentration evolution for battery simulations.
 */
class CnCH : public Concentrations {
public:
    // Constructor
    CnCH(Initialize_Geometry &geo, Domain_Parameters &para);

    // Initialize with initial concentration and psi
    void Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx);

    void TimeStep(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);

    // mfem::ParGridFunction Mub, Mob, Rxc;
    
private:
    // Geometry and domain references
    Initialize_Geometry &geometry;
    Domain_Parameters &domain_parameters;
    std::shared_ptr<mfem::ParFiniteElementSpace> fespace;

    // State vectors
    mfem::ParGridFunction Mub, Mob, Rxc;
    mfem::HypreParVector phV0, Lp1, Lp2, RHS, MuV;

    std::unique_ptr<mfem::ParBilinearForm> M_init; 
    std::shared_ptr<mfem::HypreParMatrix> MmatCH;

    // Solver and preconditioner
    std::shared_ptr<mfem::CGSolver> MCH_solver;
    mfem::HypreSmoother MCH_prec;

    std::unique_ptr<mfem::ParBilinearForm> Grad_MForm;
    std::shared_ptr<mfem::HypreParMatrix> Grad_MM; 

    std::unique_ptr<mfem::ParLinearForm> B_init; 

    std::unique_ptr<mfem::ParBilinearForm> Grad_EForm;
    std::shared_ptr<mfem::HypreParMatrix> Grad_EM;

    mfem::ParLinearForm Fct; ///< Linear form for free energy calculations
    mfem::HypreParVector Fcb; ///< Vector for storing free energy contributions

    mfem::GridFunctionCoefficient cDp;
    mfem::GridFunctionCoefficient cAp;

    std::shared_ptr<mfem::HypreParVector> CpV0; ///< Initial particle concentration values
    mfem::HypreParVector PsVc; ///< Vector for storing true degrees of freedom in the solid region
    std::shared_ptr<mfem::HypreParVector> RHCp; ///< Right-hand-side vector at the current time step
    std::shared_ptr<mfem::HypreParVector> CpVn; ///< Particle concentration values at the next time step

    // Interpolation tables
    mfem::Vector Ticks = mfem::Vector(101);
    mfem::Vector chmPot = mfem::Vector(101);
    mfem::Vector Mobility = mfem::Vector(101);
    mfem::Vector OCV = mfem::Vector(101);
    mfem::Vector i0 = mfem::Vector(101);

    // Interpolates mobility using table values
    double GetTableValues(double cn, const mfem::Vector &ticks, const mfem::Vector &data);


};

#endif // CNCH_HPP
