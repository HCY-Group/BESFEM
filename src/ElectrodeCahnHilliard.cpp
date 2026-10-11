#include "../include/ElectrodeCahnHilliard.hpp"
#include "../include/Constants.hpp"
#include "../include/MaterialProperties.hpp"
#include "../include/SimulationConfig.hpp"
#include "mfem.hpp"

ElectrodeCahnHilliard::ElectrodeCahnHilliard(Initialize_Geometry &geo, Domain_Parameters &para, sim::MaterialType mat, const SimulationConfig &cfg)
     : ConcentrationBase(geo, para, mat, cfg), cfg(cfg), pmesh(geo.parallelMesh.get()), combine_particle_groups(geo.combine_particle_groups),
      Mub(fespace.get()), Mob(fespace.get()), RxA(fespace.get()),
      Lp1(fespace.get()), Lp2(fespace.get()), MuV(fespace.get()),
      PsVc(fespace.get()), CpV0(fespace.get()), RHCp(fespace.get()), CpVn(fespace.get()),
      Fct(fespace.get()), Fcb(fespace.get()), cDp(&Mob), cAp(&RxA), MCH_solver(MPI_COMM_WORLD)

{}

void ElectrodeCahnHilliard::SetupField(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx, double gtPsx)
{
    concentration_updates = 0;
    utils.SetInitialValue(Cn, initial_value); 
    utils.CalculateLithiation(Cn, psx, gtPsx); 
    Xfr = utils.GetLithiation();

    mfem::GridFunctionCoefficient coef(&psx); 
    fem.InitializeMassMatrix(coef, M_init); 
    fem.FormSystemMatrix(M_init, boundary_dofs, MmatCH); 
    
    MCH_prec.SetType(mfem::HypreSmoother::Jacobi); 
    fem.SolverConditions(MmatCH, MCH_solver, MCH_prec); 
    // Solve the small timestep increment accurately, rather than tolerating
    // errors relative to the much larger existing lithium inventory.
    MCH_solver.SetRelTol(1.0e-12);

    // MCH_solver.iterative_mode = true;
    fem.InitializeStiffnessMatrix(cDp, Grad_MForm); 

    mfem::ConstantCoefficient varE(cfg.gc/pow(cfg.dh, pmesh->Dimension())); 
    fem.InitializeStiffnessMatrix(varE, Grad_EForm); 

    fem.InitializeForceTerm(cAp, B_init);
    Fct = *B_init;

    fem.FormLinearSystem(Grad_EForm, boundary_dofs, Cn, Fct, Grad_EM, X1v, Fcb);
    fem.FormLinearSystem(Grad_MForm, boundary_dofs, Mub, Fct, Grad_MM, X1v, Fcb); 

    psx.GetTrueDofs(PsVc); 
}

void ElectrodeCahnHilliard::UpdateConcentration(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx,
                            double gtPsx, mfem::ParGridFunction &weight_elec, const std::vector<ConcentrationBase::PairCoupling> &pair_terms)
{      
    const double rho = MaterialProperties::SiteDensity(material);
    utils.InitializeReaction(Rx, RxA, (1.0/rho)); 

    if (!combine_particle_groups)
    {
        RxA *= weight_elec;

        for (const auto &pair : pair_terms)
        {
            utils.ComputePairFlux(*pair.sum_part, *pair.weight, *pair.grad_psi, *pair.mu_self, *pair.mu_nbr, rho);
            RxA += *pair.sum_part;
        }
    }

    cAp.SetGridFunction(&RxA); 

    fem.Update(B_init); 
    Fct = *B_init; 

    // Tabulate the chemical potential and mobility values
    for (int i = 0; i < Cn.Size(); i++)
    {
        // double val = Cn(i);
        double val = std::clamp(Cn(i), 1.0e-6, 0.96);

        double mu = MaterialProperties::ChemicalPotential(material, val);

        if (std::abs(mu) > 1.0e4)
        {
            mu = ((mu / (-1*Constants::Frd)) - 3.4) * -1;
        }

        Mub(i) = mu;
        Mob(i) = psx(i) * MaterialProperties::Mobility(material, val);
    }

    Cn.GetTrueDofs(CpV0);

    // Lap phi
    Grad_EM.Mult(CpV0, Lp1); 
    Mub.GetTrueDofs(MuV); 
    MuV += Lp1; 

    // Update stiffness form with new mobility
    cDp.SetGridFunction(&Mob); 
    fem.Update(Grad_MForm);
    fem.FormLinearSystem(Grad_MForm, boundary_dofs, Mub, Fct, Grad_MM, X1v, Fcb); 

    Grad_MM.Mult(MuV, Lp2); 
    Lp2.Neg(); 
    Lp2 *= cfg.dt; 

    // Add the reaction term to the right-hand side vector
    Fcb *= cfg.dt;
    Lp2 += Fcb; 

    // M * delta_c = dt * (reaction - transport). Keeping the baseline out
    // of the solve prevents CG error from acting as a per-step lithium sink.
    RHCp = Lp2;
    MCH_solver.Mult(RHCp, CpVn); 
    MFEM_VERIFY(MCH_solver.GetConverged(), "CH concentration increment solve failed.");
    CpVn += CpV0;

    const bool diagnose = concentration_updates++ % std::max(1, cfg.save_freq) == 0;
    double balance[3] = {}; // raw inventory increment, RHS, clipping change
    if (diagnose)
    {
        mfem::Vector delta(CpVn), inventory(CpVn.Size());
        delta -= CpV0;
        MmatCH.Mult(delta, inventory);
        balance[0] = inventory.Sum();
        balance[1] = Lp2.Sum();
    }
    mfem::Vector raw(CpVn);
    long long clips[2] = {};

    // Ensure that the concentration values are within the valid range
    for (int i = 0; i < CpV0.Size(); i++) {
        // The phase mask has a positive floor in Domain_Parameters. Its tail
        // therefore has storage and must retain its update for conservation.

        if (CpVn(i) < 0.0){
            ++clips[0];
            (CpVn)(i) = 1.0e-6;}

        if (CpVn(i) > 0.96){
            ++clips[1];
            (CpVn)(i) = 0.96;}
    }

    if (diagnose)
    {
        mfem::Vector correction(CpVn), inventory(CpVn.Size());
        correction -= raw;
        MmatCH.Mult(correction, inventory);
        balance[2] = inventory.Sum();
        MPI_Allreduce(MPI_IN_PLACE, balance, 3, MPI_DOUBLE, MPI_SUM, fespace->GetComm());
        MPI_Allreduce(MPI_IN_PLACE, clips, 2, MPI_LONG_LONG_INT, MPI_SUM, fespace->GetComm());
        if (mfem::Mpi::WorldRank() == 0)
        {
            std::cout << "[CH balance] update=" << concentration_updates
                      << " raw_delta=" << balance[0] << " expected_delta=" << balance[1]
                      << " clipping_delta=" << balance[2]
                      << " lower_clips=" << clips[0] << " upper_clips=" << clips[1] << std::endl;
        }
    }

    // recover the GridFunction from the HypreParVector
    Cn.Distribute(CpVn); 

    utils.CalculateLithiation(Cn, psx, gtPsx); 
    Xfr = utils.GetLithiation();

    // Rx = RxA; 
}



