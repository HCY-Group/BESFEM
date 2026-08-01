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
    utils.SetInitialValue(Cn, initial_value); 
    utils.CalculateLithiation(Cn, psx, gtPsx); 
    Xfr = utils.GetLithiation();


    mfem::GridFunctionCoefficient coef(&psx); 
    fem.InitializeMassMatrix(coef, M_init); 
    fem.FormSystemMatrix(M_init, boundary_dofs, MmatCH); 
    
    // MCH_solver.iterative_mode = false; 
    MCH_prec.SetType(mfem::HypreSmoother::Jacobi); 
    fem.SolverConditions(MmatCH, MCH_solver, MCH_prec); 

    MCH_solver.iterative_mode = true;
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

    if (!combine_particle_groups){
        
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

        // Mub(i) = mu;


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

    // Update the right-hand side vector for the Cahn-Hilliard equation
    MmatCH.Mult(CpV0, RHCp); 
    RHCp += Lp2; 

    MCH_solver.Mult(RHCp, CpVn); 

    // Ensure that the concentration values are within the valid range
    for (int i = 0; i < CpV0.Size(); i++) {
        if (PsVc(i) < 1.0e-5) {
            (CpVn)(i) = CpV0(i);; 
        }

        if (CpVn(i) < 0.0){
            (CpVn)(i) = 1.0e-6;}

        if (CpVn(i) > 0.96){
            (CpVn)(i) = 0.96;}
    }

    // recover the GridFunction from the HypreParVector
    Cn.Distribute(CpVn); 

    utils.CalculateLithiation(Cn, psx, gtPsx); 
    Xfr = utils.GetLithiation();

    // Rx = RxA; 
}



