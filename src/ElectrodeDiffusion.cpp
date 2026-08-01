#include "../include/ElectrodeDiffusion.hpp"
#include "../include/Constants.hpp"
#include "../include/MaterialProperties.hpp"
#include "../include/SimulationConfig.hpp"
#include "../include/SimTypes.hpp"
#include "mfem.hpp"

ElectrodeDiffusion::ElectrodeDiffusion(Initialize_Geometry &geo, Domain_Parameters &para, sim::MaterialType mat, const SimulationConfig &cfg)
    : ConcentrationBase(geo, para, mat, cfg),  cfg(cfg), combine_particle_groups(geo.combine_particle_groups), Dp(fespace.get()), Mp_solver(MPI_COMM_WORLD), cDp(&Dp), Fct(fespace.get()), PsVc(fespace.get()), CpV0(fespace.get()), RHCp(fespace.get()), 
    CpVn(fespace.get()), Rxn(fespace.get()), cAp(&Rxn)
{}

void ElectrodeDiffusion::SetupField(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx, double gtPsx)
{
    boundary_dofs.DeleteAll();

    utils.SetInitialValue(Cn, initial_value);
    utils.CalculateLithiation(Cn, psx, gtPsx);
    Xfr = utils.GetLithiation();

    if (mfem::Mpi::WorldRank() == 0)
    {
        std::cout << "Initial Lithiation Fraction: " << Xfr << std::endl;
    }

    mfem::GridFunctionCoefficient coef(&psx);
    fem.InitializeMassMatrix(coef, Mt);
    fem.FormSystemMatrix(Mt, boundary_dofs, Mmatp);

    Mp_solver.iterative_mode = false;
    Mp_prec.SetType(mfem::HypreSmoother::Jacobi); 
    fem.SolverConditions(Mmatp, Mp_solver, Mp_prec);
    
    fem.InitializeStiffnessMatrix(cDp, Kc2);
    psx.GetTrueDofs(PsVc);

}

void ElectrodeDiffusion::UpdateConcentration(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx,
                        double gtPsx, mfem::ParGridFunction &weight_elec, const std::vector<ConcentrationBase::PairCoupling> &pair_terms)

{
    const double rho = MaterialProperties::SiteDensity(material);
    utils.InitializeReaction(Rx, Rxn, (1.0/rho));

    if (!combine_particle_groups){
        Rxn *= weight_elec;
        for (const auto &pair : pair_terms)
        {
            utils.ComputePairFlux(*pair.sum_part, *pair.weight, *pair.grad_psi, *pair.mu_self, *pair.mu_nbr, rho);
            Rxn += *pair.sum_part;
        }
    }

    cAp.SetGridFunction(&Rxn);

    fem.InitializeForceTerm(cAp, Bc2);
    fem.Update(Bc2);
    Fct = *Bc2;

    for (int vi = 0; vi < nV; vi++){
        const double cn_val = Cn(vi);
        Dp(vi) = psx(vi) * MaterialProperties::Diffusivity(material, cn_val);
        if (Dp(vi) > 4.6e-10){Dp(vi) = 4.6e-10;}
    }
    cDp.SetGridFunction(&Dp);

    fem.Update(Kc2);
    fem.FormLinearSystem(Kc2, boundary_dofs, Cn, Fct, Kmatp, X1v, Fcb);
    Fcb *= cfg.dt;

    Tmatp.reset(Add(1.0, Mmatp, -1.0* cfg.dt, Kmatp));

    Cn.GetTrueDofs(CpV0);
    Tmatp->Mult(CpV0, RHCp);
    RHCp += Fcb;

    Mp_solver.Mult(RHCp, CpVn);

    psx.GetTrueDofs(PsVc);

    for (int p = 0; p < CpV0.Size(); p++){
        if (PsVc(p) < 1.0e-5){ 
            (CpVn)(p) = CpV0(p);} 

        if (CpVn(p) < 0.0){
            (CpVn)(p) = 1.0e-6;}
    }

    Cn.Distribute(CpVn);

    utils.CalculateLithiation(Cn, psx, gtPsx);
    Xfr = utils.GetLithiation();

    Rx = Rxn;

}