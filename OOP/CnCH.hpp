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

private:
    // Geometry and domain references
    Initialize_Geometry &geometry;
    Domain_Parameters &domain_parameters;
    std::shared_ptr<mfem::ParFiniteElementSpace> fespace;

    // State vectors
    mfem::ParGridFunction Mub, Mob;
    mfem::HypreParVector phV0, Lp1, Lp2, RHS, MuV;

    mfem::ParGridFunction AvP; ///< Grid function for active particle surface area

    // Matrices
    std::shared_ptr<mfem::HypreParMatrix> M_phi;
    std::shared_ptr<mfem::HypreParMatrix> K_phi, K_mu;

    // Solver and preconditioner
    std::shared_ptr<mfem::CGSolver> MCH_solver;
    mfem::HypreSmoother MCH_prec;

    std::shared_ptr<mfem::HypreParVector> CpV0; ///< Initial particle concentration values
    mfem::HypreParVector PsVc; ///< Vector for storing true degrees of freedom in the solid region
    std::shared_ptr<mfem::HypreParVector> RHCp; ///< Right-hand-side vector at the current time step
    std::shared_ptr<mfem::HypreParVector> CpVn; ///< Particle concentration values at the next time step


    // Interpolation tables
    mfem::Vector X_101 = mfem::Vector(101);
    mfem::Vector dF_101 = mfem::Vector(101);
    mfem::Vector Mb5_101 = mfem::Vector(101);

    /// Interpolates dF using table values
    double Calculate_dF(double ph, const mfem::Vector &X_101, const mfem::Vector &dF_101);

    /// Interpolates mobility using table values
    double Calculate_Mobility(double cn, const mfem::Vector &X_101, const mfem::Vector &mob_101);


};

#endif // CNCH_HPP
