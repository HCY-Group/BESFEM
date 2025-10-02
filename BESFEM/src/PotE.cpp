/**
 * @file PotE.cpp
 * @brief Implementation of the potential class for electrolyte potential simulations.
 */

#include "../include/PotE.hpp"
#include "../inputs/Constants.hpp"
#include "mfem.hpp"
#include "../include/CnE.hpp"
#include <optional>


PotE::PotE(Initialize_Geometry &geo, Domain_Parameters &para)
    : Potentials(geo,para), geometry(geo), domain_parameters(para), fespace(geo.parfespace), dbc_bdr(geo.dbc_bdr), gtPse(para.gtPse), 
    kpl(fespace.get()), RpE(fespace.get()), Dmp(fespace.get()), pE0(fespace.get())
    
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

// Make center = BvE, like the Fortran gauge fix.
static inline void GaugeToCenterValue(mfem::ParGridFunction &phx, double BvE)
{
    auto *pfes = phx.ParFESpace();
    auto &mesh = *pfes->GetParMesh();  // ParMesh is fine; Mesh also works
    const int dim = mesh.Dimension();

    // 1) Global geometric center (bounding box mid-point)
    mfem::Vector xmin(dim), xmax(dim), xcen(dim);
    mesh.GetBoundingBox(xmin, xmax);
    for (int d = 0; d < dim; ++d) xcen(d) = 0.5 * (xmin(d) + xmax(d));

    // 2) Find locally closest element by centroid distance^2
    double best_loc_d2 = std::numeric_limits<double>::infinity();
    int    best_loc_el = -1;
    mfem::Vector c(dim);
    for (int el = 0; el < mesh.GetNE(); ++el)
    {
        mfem::ElementTransformation *T = mesh.GetElementTransformation(el);
        // Centroid in ref. coords is (0,0[,0]); map to physical
        mfem::IntegrationPoint ip; ip.Set3(0.0, 0.0, 0.0); // works for 2D/3D; z ignored in 2D
        T->SetIntPoint(&ip);
        mfem::Vector x(dim);
        T->Transform(ip, x);

        double d2 = 0.0;
        for (int d = 0; d < dim; ++d) { double dd = x(d) - xcen(d); d2 += dd*dd; }
        if (d2 < best_loc_d2) { best_loc_d2 = d2; best_loc_el = el; c = x; }
    }

    // 3) Pick the globally closest element with MINLOC
    struct { double val; int rank; } in{best_loc_d2, pfes->GetMyRank()}, out{0.0, 0};
    MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, pfes->GetComm());

    // 4) The winner rank evaluates φ_e at its best element's centroid
    double ph_center = 0.0;
    if (pfes->GetMyRank() == out.rank)
    {
        mfem::ElementTransformation *T = mesh.GetElementTransformation(best_loc_el);
        mfem::IntegrationPoint ip; ip.Set3(0.0, 0.0, 0.0);
        T->SetIntPoint(&ip);
        ph_center = phx.GetValue(*T, ip);
    }

    // 5) Broadcast φ(center) and compute DvPe = φ(center) - BvE
    MPI_Bcast(&ph_center, 1, MPI_DOUBLE, out.rank, pfes->GetComm());
    const double DvPe = ph_center - BvE;

    // 6) Apply the constant shift φ := φ - DvPe   (now φ(center) == BvE)
    phx += -DvPe;
}


void PotE::Initialize(mfem::ParGridFunction &ph, double initial_value, mfem::ParGridFunction &psx)
{
    BvE = initial_value; // Set the boundary value.
    Potentials::SetInitialPotentials(ph, BvE); // Initialize potentials

    SolverSteps::InitializeStiffnessMatrix(cKe, Kl2); // Initialize the stiffness matrix

    mfem::ConstantCoefficient dbc_potE_Coef(BvE); // Coefficient for Dirichlet boundary conditions HALF
    ph.ProjectBdrCoefficient(dbc_potE_Coef, dbc_bdr); // Apply Dirichlet boundary conditions HALF

    fespace->GetEssentialTrueDofs(dbc_bdr, ess_tdof_list_potE); // Get essential true degrees of freedom for Dirichlet boundary conditions HALF

    SolverSteps::FormLinearSystem(Kl2, ess_tdof_list_potE, ph, B1t, Kml, X1v, B1v); // Assemble the linear system HALF
    // SolverSteps::FormLinearSystem(Kl2, boundary_dofs, ph, B1t, Kml, X1v, B1v); // Assemble the linear system FULL

    Mpe = std::make_unique<mfem::HypreBoomerAMG>(Kml);  // builds hierarchy once
    Mpe->SetPrintLevel(0);
    SolverSteps::SolverConditions(Kml, cgPE_solver, *Mpe); // Set up the solver conditions

    SolverSteps::InitializeForceTerm(cRe, Bl2); // Initialize the force term
    SolverSteps::Update(Bl2); // Update the force term
    Flt = *Bl2; // Move the force term

    SolverSteps::FormLinearSystem(Kl2, ess_tdof_list_potE, ph, Flt, Kml, X1v, Flb); // Assemble the force term system HALF
    // SolverSteps::FormLinearSystem(Kl2, boundary_dofs, ph, Flt, Kml, X1v, Flb); // Assemble the force term system FULL


    SolverSteps::InitializeStiffnessMatrix(cDm, Kl1); // Initialize the diffusivity matrix
    SolverSteps::FormLinearSystem(Kl1, boundary_dofs, ph, B1t, Kdm, X1v, B1v); // Assemble the diffusivity matrix system
 

}

void PotE::TimeStep(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &potential, mfem::HypreParVector &CeVn)
{
    ElectrolyteConductivity(Cn, psx); // Update conductivity and diffusivity

    SolverSteps::Update(Kl1); // Update the diffusivity matrix
    SolverSteps::FormLinearSystem(Kl1, boundary_dofs, potential, B1t, Kdm, X1v, B1v); // Assemble the diffusivity matrix system

    Cn.GetTrueDofs(CeVn); // Get the true degrees of freedom for concentration
    Kdm.Mult(CeVn, LpCe); 

    // std::cout << "LpCe in TimeStep: " << LpCe.Sum() << std::endl;

    SolverSteps::Update(Kl2); // Update the conductivity matrix

    mfem::ConstantCoefficient dbc_potE_Coef(BvE); // Coefficient for Dirichlet boundary conditions HALF
    potential.ProjectBdrCoefficient(dbc_potE_Coef, dbc_bdr); // Apply Dirichlet boundary conditions HALF

    SolverSteps::FormLinearSystem(Kl2, ess_tdof_list_potE, potential, B1t, Kml, X1v, B1v); // Assemble the conductivity matrix system HALF
    // SolverSteps::FormLinearSystem(Kl2, boundary_dofs, potential, B1t, Kml, X1v, B1v); // Assemble the conductivity matrix system FULL


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

// void PotE::Advance(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror)
// {
//     RpE = Rx1;
//     RpE += Rx2;

//     RpE.Neg(); // check and see if this should be negative
    
//     Bl2->Assemble();
//     Flt = *Bl2;
    
//     SolverSteps::FormLinearSystem(Kl2, boundary_dofs, phx, Flt, Kml, X1v, Flb); // Assemble the force term system

//     RHSl = Flb;
//     RHSl += LpCe;

//     pE0 = phx; // Store the current potential field
//     pE0.GetTrueDofs(Xe0); // Extract degrees of freedom

//     phx.Save("phE_gf_before");

//     cgPE_solver.Mult(RHSl, Xe0); // Solve for the error term

//     std::cout << Xe0.Sum() << " Xe0 after solver" << std::endl;

//     phx.Distribute(Xe0); // Distribute the updated values
//     GaugeToCenterValue(phx, BvE);


//     phx.Save("phE_gf_after");

//     Potentials::ComputeGlobalError(pE0, phx, psx, gerror, gtPse); // Compute global error
// }

void PotE::Advance(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2,
                   mfem::ParGridFunction &phx, mfem::ParGridFunction &psx,
                   double &gerror)
{
    RpE = Rx1;
    RpE += Rx2;
    RpE.Neg(); // keeps your original sign choice

    Bl2->Assemble();
    Flt = *Bl2;

    // Assemble A, rhs parts as you already do
    SolverSteps::FormLinearSystem(Kl2, boundary_dofs, phx, Flt, Kml, X1v, Flb);

    RHSl = Flb;
    RHSl += LpCe;

    // ---- Enforce Neumann compatibility: (rhs, 1) = 0 ----
    mfem::ParLinearForm one_lf(fespace.get());
    mfem::ConstantCoefficient one_coef(1.0);
    one_lf.AddDomainIntegrator(new mfem::DomainLFIntegrator(one_coef));
    one_lf.Assemble();
    // true-dof vector corresponding to constant-1 test
    mfem::HypreParVector *w_const = one_lf.ParallelAssemble();

    double num = mfem::InnerProduct(RHSl, *w_const);
    double den = mfem::InnerProduct(*w_const, *w_const);
    if (den > 0.0) { RHSl.Add(-num/den, *w_const); }

    delete w_const;
    // -----------------------------------------------------

    // save current field as initial guess
    pE0 = phx;
    pE0.GetTrueDofs(Xe0);
    phx.Save("phE_gf_before");

    // --------- STRICTLY-SPD PRECONDITIONER FOR CG ----------
    mfem::DSmoother J(0);      // 0 = Jacobi (SPD)
    J.SetOperator(Kml);        // attach the matrix

    cgPE_solver.SetOperator(Kml);
    cgPE_solver.SetPreconditioner(J);   // <-- use DSmoother, not AMG here
    cgPE_solver.iterative_mode = true; // use Xe0 as initial guess
    cgPE_solver.SetMaxIter(500);
    cgPE_solver.SetRelTol(1e-8);
    cgPE_solver.SetAbsTol(0.0);
    cgPE_solver.SetPrintLevel(1);
    // -------------------------------------------------------

    cgPE_solver.Mult(RHSl, Xe0);

    std::cout << Xe0.Sum() << " Xe0 after solver" << std::endl;

    phx.Distribute(Xe0);
    GaugeToCenterValue(phx, BvE);          // keep your gauge
    phx.Save("phE_gf_after");

    Potentials::ComputeGlobalError(pE0, phx, psx, gerror, gtPse);


}



void PotE::ElectrolyteConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx) {
    for (int vi = 0; vi < nV; vi++){
        dffe = exp(-7.02 - 830 * Cn(vi) + 50000 * Cn(vi) * Cn(vi)); // Compute diffusivity factor
        Dmp(vi) = psx(vi) * tc1 * Constants::D0 * dffe;
        kpl(vi) = psx(vi) * tc2 * Constants::D0 * dffe * Cn(vi);
    }
}
