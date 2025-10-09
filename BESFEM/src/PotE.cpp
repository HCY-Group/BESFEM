/**
 * @file PotE.cpp
 * @brief Implementation of the potential class for electrolyte potential simulations.
 */

#include "../include/PotE.hpp"
#include "../inputs/Constants.hpp"
#include "mfem.hpp"
#include "../include/CnE.hpp"
#include <optional>

// PickElectrolyteAnchorTDof.hpp
// #pragma once
#include "mfem.hpp"
#include <limits>

// /// Returns the *local* true-dof index (on the owning rank) of a single electrolyte
// /// anchor tdof to pin. Returns -1 on ranks that do not own it.
// /// - fes    : ParFiniteElementSpace for phi_E  (NON-const! MFEM ctor needs non-const)
// /// - pse    : electrolyte mask (~1 in electrolyte). Pass nullptr to skip mask.
// /// - thresh : choose tdof with pse > thresh (e.g., 0.8–0.95 to avoid interface)
// static inline int PickElectrolyteAnchorTDof(mfem::ParFiniteElementSpace &fes,
//                                             const mfem::ParGridFunction *pse,
//                                             double thresh = 0.5)
// {
//     MPI_Comm comm = fes.GetComm();
//     int rank = 0, nprocs = 1;
//     MPI_Comm_rank(comm, &rank);
//     MPI_Comm_size(comm, &nprocs);

//     // Global partitioning of true DOFs: offsets[0..nprocs] with global range [g0,g1)
//     const HYPRE_BigInt *offsets = fes.GetTrueDofOffsets(); // length nprocs+1
//     const HYPRE_BigInt g0 = offsets[rank];
//     const HYPRE_BigInt g1 = offsets[rank+1];
//     const HYPRE_BigInt Nglob = offsets[nprocs];


//     // Local true DOF count (MFEM returns int here)
//     const int nloc = fes.TrueVSize();

//     int local_idx = -1;
//     if (pse && nloc > 0) {
//         mfem::HypreParVector pse_tdof(&fes);
//         pse->GetTrueDofs(pse_tdof);
//         const double *d = pse_tdof.GetData();
//         for (int i = 0; i < nloc; i++) {
//             if (d[i] > thresh) { local_idx = i; break; }
//         }
//     }

//     // B) turn local candidate into a comparable GLOBAL index
//     const HYPRE_BigInt BIG = std::numeric_limits<HYPRE_BigInt>::max();
//     HYPRE_BigInt my_global = (local_idx >= 0) ? (g0 + (HYPRE_BigInt)local_idx) : BIG;

//     // C) choose a single global winner (lowest global tdof where pse>thresh)
//     HYPRE_BigInt winner_global = BIG;
//     MPI_Allreduce(&my_global, &winner_global, 1, HYPRE_MPI_BIG_INT, MPI_MIN, comm);

//     // D) robust fallback: if no candidate anywhere, pick the first rank with any dofs
//     if (winner_global == BIG) {
//         int my_rank_if_has = (nloc > 0) ? rank : INT_MAX;
//         int owner_rank = INT_MAX;
//         MPI_Allreduce(&my_rank_if_has, &owner_rank, 1, MPI_INT, MPI_MIN, comm);
//         if (owner_rank == INT_MAX || Nglob == 0) { return -1; } // degenerate case
//         return (rank == owner_rank) ? 0 : -1; // pin local tdof 0 on the first nonempty rank
//     } 

//     // E) find the owner rank of winner_global using offsets
//     int owner_rank = -1;
//     // (binary search would be fine, but linear is okay for small nprocs)
//     for (int r = 0; r < nprocs; ++r) {
//         if (winner_global >= offsets[r] && winner_global < offsets[r+1]) { owner_rank = r; break; }
//     }

//     // F) make everyone agree on the owner: use MAX (only owner has nonnegative rank)
//     int i_own = (owner_rank == rank) ? rank : -1;
//     int agreed_owner = -1;
//     MPI_Allreduce(&i_own, &agreed_owner, 1, MPI_INT, MPI_MAX, comm);

//     // If somehow ambiguous, fall back to first rank with dofs (shouldn’t happen)
//     if (agreed_owner < 0) {
//         int my_rank_if_has = (nloc > 0) ? rank : INT_MAX;
//         MPI_Allreduce(&my_rank_if_has, &agreed_owner, 1, MPI_INT, MPI_MIN, comm);
//         if (agreed_owner == INT_MAX) return -1;
//         return (rank == agreed_owner) ? 0 : -1;
//     }

//     if (rank == agreed_owner) {
//         // local index = winner_global - g0 (for this owner)
//         return (int)(winner_global - offsets[rank]);
//     }
//     return -1;
// }

// #include <climits>

static inline int PickAnchor(mfem::ParFiniteElementSpace &fes, HYPRE_BigInt g)
{
    // Global true-dof ownership ranges
    const HYPRE_BigInt *offs = fes.GetTrueDofOffsets(); // size = nprocs+1 - perhaps work for parallel?
    int rank = 0, nprocs = 0;
    MPI_Comm_rank(fes.GetComm(), &rank); 
    MPI_Comm_size(fes.GetComm(), &nprocs);

    // Defensive: space must be finalized (have true dofs)
    const HYPRE_BigInt n_global_tdofs = offs ? offs[nprocs] : 0;
    if (n_global_tdofs == 0) {
        if (rank == 0) {
            std::cerr << "[Anchor] ERROR: ParFiniteElementSpace has zero true DOFs.\n"
                         "         Construct/Update the space before picking an anchor.\n";
        }
        return -1;
    }

    // check requested g is in range
    if (g < 0 || g >= offs[nprocs]) {
        if (rank == 0) {
            std::cerr << "[Anchor] Requested global tdof " << g
                      << " is out of range [0, " << (offs[nprocs]-1) << "].\n";
        }
        return -1;
    }

    // Does this rank own g? Need to ensure exactly one rank does this so that the DOF can be set
    if (g >= offs[rank] && g < offs[rank+1]) {
        const HYPRE_BigInt local = g - offs[rank];
        std::cout << "[Rank " << rank << "] Anchoring GLOBAL tdof " << g
                  << " (local index " << local << ")\n";
        return static_cast<int>(local);
    }

    // All other ranks: no anchor here
    return -1;
}

PotE::PotE(Initialize_Geometry &geo, Domain_Parameters &para)
    : Potentials(geo,para), geometry(geo), domain_parameters(para), fespace(geo.parfespace), dbc_bdr(geo.dbc_bdr), gtPse(para.gtPse), 
    kpl(fespace.get()), RpE(fespace.get()), Dmp(fespace.get()), pE0(fespace.get()), phE_bc(fespace.get())
    
    {
    cgPE_solver = mfem::CGSolver(MPI_COMM_WORLD);
    
    B1t = mfem::ParLinearForm(fespace.get());
    X1v = mfem::HypreParVector(fespace.get());
    B1v = mfem::HypreParVector(fespace.get());
    Flb = mfem::HypreParVector(fespace.get());
    LpCe = mfem::HypreParVector(fespace.get()); // Initialize the vector for concentration degrees of freedom
    RpE = mfem::ParGridFunction(fespace.get());
    Dmp = mfem::ParGridFunction(fespace.get()); // Initialize diffusivity field

    Bl2 = std::make_unique<mfem::ParLinearForm>(fespace.get());

    kpl = mfem::ParGridFunction(fespace.get()); // Initialize conductivity field
    cKe = mfem::GridFunctionCoefficient(&kpl); // Coefficient for conductivity field
    cRe = mfem::GridFunctionCoefficient(&RpE);
    Flt = mfem::ParLinearForm(fespace.get());
    cDm = mfem::GridFunctionCoefficient(&Dmp); // Coefficient for diffusivity field

    Kl1 = std::make_unique<mfem::ParBilinearForm>(fespace.get()); // Initialize the bilinear form for conductivity
    Kl2 = std::make_unique<mfem::ParBilinearForm>(fespace.get()); // Initialize the bilinear form for conductivity

    pE0 = mfem::ParGridFunction(fespace.get()); // Initialize the potential grid function
    Xe0 = mfem::HypreParVector(fespace.get()); // Initialize the solution vector for potential
    RHSl = mfem::HypreParVector(fespace.get()); // Initialize the right-hand side vector for potential

    }

#include <mpi.h>

void PotE::Initialize(mfem::ParGridFunction &ph, double initial_value, mfem::ParGridFunction &psx)
{

    // const HYPRE_BigInt anchor_global = 4;

    
    // // --- pick anchor ---
    // if (!anchor_ready) {
    //     const int a = PickAnchor(*fespace, anchor_global); // pick global tdof to anchor
    //     if (a >= 0) { ess_tdof_potE.SetSize(1); ess_tdof_potE[0] = a; } // if tdof on this rank, set ess_tdof_potE to it
    //     else        { ess_tdof_potE.SetSize(0); } // otherwise, no anchor on this rank
    //     anchor_ready = true; 
    // }

    const int anchor = 2000;

    if (!anchor_set) {
        ess_tdof_potE.SetSize(1);
        ess_tdof_potE[0] = anchor;  
        anchor_set = true; 
    }


    
    BvE = initial_value; // Set the boundary value.
    Potentials::SetInitialPotentials(ph, BvE); // Initialize potentials

    // --- BC vector for gauge ---
    phE_bc = ph;                       // start from current iterate
    if (ess_tdof_potE.Size() == 1) {
        mfem::Vector td; 
        phE_bc.GetTrueDofs(td); // get true dofs
        td(ess_tdof_potE[0]) = Constants::init_BvE;     // set anchor value
        phE_bc.SetFromTrueDofs(td); // set grid function from modified true dofs
    }
    
    SolverSteps::InitializeStiffnessMatrix(cKe, Kl2); // Initialize the stiffness matrix

    // mfem::ConstantCoefficient dbc_potE_Coef(BvE); // Coefficient for Dirichlet boundary conditions HALF
    // ph.ProjectBdrCoefficient(dbc_potE_Coef, dbc_bdr); // Apply Dirichlet boundary conditions HALF

    // fespace->GetEssentialTrueDofs(dbc_bdr, ess_tdof_list_potE); // Get essential true degrees of freedom for Dirichlet boundary conditions HALF

    // std::cout << "Number of essential true dofs for PotE: " << ess_tdof_list_potE.Size() << std::endl;
    // std::cout << "boundary_dofs size for PotE: " << boundary_dofs.Size() << std::endl;

    // SolverSteps::FormLinearSystem(Kl2, ess_tdof_list_potE, ph, B1t, Kml, X1v, B1v); // Assemble the linear system HALF
    SolverSteps::FormLinearSystem(Kl2, ess_tdof_potE, phE_bc, B1t, Kml, X1v, B1v); // Assemble the linear system FULL

    std::cout << "Linear system for PotE assembled." << std::endl;

    Mpe = std::make_unique<mfem::HypreBoomerAMG>(Kml);  // builds hierarchy once
    Mpe->SetPrintLevel(0);
    SolverSteps::SolverConditions(Kml, cgPE_solver, *Mpe); // Set up the solver conditions

    std::cout << "Solver conditions for PotE set." << std::endl;

    SolverSteps::InitializeForceTerm(cRe, Bl2); // Initialize the force term
    SolverSteps::Update(Bl2); // Update the force term
    Flt = *Bl2; // Move the force term

    std::cout << "Force term for PotE initialized." << std::endl;

    // SolverSteps::FormLinearSystem(Kl2, ess_tdof_list_potE, ph, Flt, Kml, X1v, Flb); // Assemble the force term system HALF
    SolverSteps::FormLinearSystem(Kl2, ess_tdof_potE, phE_bc, Flt, Kml, X1v, Flb); // Assemble the force term system FULL

    std::cout << "Force term system for PotE assembled." << std::endl;

    SolverSteps::InitializeStiffnessMatrix(cDm, Kl1); // Initialize the diffusivity matrix
    SolverSteps::FormLinearSystem(Kl1, ess_tdof_potE, phE_bc, B1t, Kdm, X1v, B1v); // Assemble the diffusivity matrix system

    std::cout << "Diffusivity matrix for PotE initialized." << std::endl;
 

}

void PotE::TimeStep(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &potential, mfem::HypreParVector &CeVn)
{
    ElectrolyteConductivity(Cn, psx); // Update conductivity and diffusivity

    SolverSteps::Update(Kl1); // Update the diffusivity matrix
    SolverSteps::FormLinearSystem(Kl1, ess_tdof_potE, phE_bc, B1t, Kdm, X1v, B1v); // Assemble the diffusivity matrix system

    Cn.GetTrueDofs(CeVn); // Get the true degrees of freedom for concentration
    Kdm.Mult(CeVn, LpCe); 

    SolverSteps::Update(Kl2); // Update the conductivity matrix

    // mfem::ConstantCoefficient dbc_potE_Coef(BvE); // Coefficient for Dirichlet boundary conditions HALF
    // potential.ProjectBdrCoefficient(dbc_potE_Coef, dbc_bdr); // Apply Dirichlet boundary conditions HALF

    // SolverSteps::FormLinearSystem(Kl2, ess_tdof_list_potE, potential, B1t, Kml, X1v, B1v); // Assemble the conductivity matrix system HALF
    SolverSteps::FormLinearSystem(Kl2, ess_tdof_potE, phE_bc, B1t, Kml, X1v, B1v); // Assemble the conductivity matrix system FULL

    Mpe->SetOperator(Kml); // Set the operator for the preconditioner
    cgPE_solver.SetPreconditioner(*Mpe);
    cgPE_solver.SetOperator(Kml); // Set the operator for the solver

}

void PotE::Advance(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror)
{
    RpE = Rx;
    RpE.Neg();

    Bl2->Assemble();
    Flt = *Bl2;

    mfem::ConstantCoefficient dbc_potE_Coef(BvE); // Coefficient for Dirichlet boundary conditions
    phx.ProjectBdrCoefficient(dbc_potE_Coef, dbc_bdr); // Apply Dirichlet boundary conditions
    
    SolverSteps::FormLinearSystem(Kl2, ess_tdof_list_potE, phx, Flt, Kml, X1v, Flb); // Assemble the force term system

    RHSl = Flb;
    RHSl += LpCe;

    pE0 = phx; // Store the current potential field
    pE0.GetTrueDofs(Xe0); // Extract degrees of freedom

    cgPE_solver.Mult(RHSl, Xe0); // Solve for the error term

    phx.Distribute(Xe0); // Distribute the updated values

    Potentials::ComputeGlobalError(pE0, phx, psx, gerror, gtPse); // Compute global error
}

void PotE::Advance(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror)
{
    RpE = Rx1;
    RpE += Rx2;

    RpE.Neg(); // check and see if this should be negative
    
    Bl2->Assemble();
    Flt = *Bl2;

    // Enforce the same gauge
    phE_bc = phx;
    if (ess_tdof_potE.Size() == 1) {
        mfem::Vector td; phE_bc.GetTrueDofs(td);
        td(ess_tdof_potE[0]) = BvE;
        phE_bc.SetFromTrueDofs(td);
    }

    // std::cout << "Number of essential true dofs for PotE: " << ess_tdof_potE.Size() << std::endl;
    SolverSteps::FormLinearSystem(Kl2, ess_tdof_potE, phE_bc, Flt, Kml, X1v, Flb);

    RHSl = Flb;
    RHSl += LpCe;

    pE0 = phx; // Store the current potential field
    pE0.GetTrueDofs(Xe0); // Extract degrees of freedom

    cgPE_solver.Mult(RHSl, Xe0); // Solve for the error term

    phx.Distribute(Xe0); // Distribute the updated values
    Potentials::ComputeGlobalError(pE0, phx, psx, gerror, gtPse); // Compute global error
}




void PotE::ElectrolyteConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx) {
    for (int vi = 0; vi < nV; vi++){
        dffe = exp(-7.02 - 830 * Cn(vi) + 50000 * Cn(vi) * Cn(vi)); // Compute diffusivity factor
        Dmp(vi) = psx(vi) * tc1 * Constants::D0 * dffe;
        kpl(vi) = psx(vi) * tc2 * Constants::D0 * dffe * Cn(vi);
    }
}
