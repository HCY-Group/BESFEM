#include "../include/ElectrolyteDiffusion.hpp"
#include "../include/Constants.hpp"
#include "../include/MaterialProperties.hpp"
#include "../include/SimulationConfig.hpp"
#include "mfem.hpp"

ElectrolyteDiffusion::ElectrolyteDiffusion(Initialize_Geometry &geo, Domain_Parameters &para, BoundaryConditions &bc, sim::CellMode mode, sim::MaterialType mat, const SimulationConfig &cfg)
    : ConcentrationBase(geo, para, mat, cfg), cfg(cfg), boundary_conditions(bc), De(fespace.get()), Rxe(fespace.get()), PeR(fespace.get()), cDe(&De), cAe(&Rxe), matCoef_R(&PeR), nbc_bdr(bc.nbc_bdr), 
    nbcCoef(0.0), Fet(fespace.get()), Me_solver(MPI_COMM_WORLD),
    CeV0(fespace.get()), RHCe(fespace.get()), CeVn(fespace.get()),
    Feb(fespace.get()), X1v(fespace.get()), PsVc(fespace.get()), mode_(mode)
    
{}

void ElectrolyteDiffusion::SetupField(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx, double gtPsx)
{
    if(mode_ == sim::CellMode::HALF){
       if (mfem::Mpi::WorldRank() == 0) {std::cout << "Initializing CnE for Half Cell" << std::endl;}

        utils.SetInitialValue(Cn, initial_value);

        mfem::GridFunctionCoefficient coef(&psx);
        fem.InitializeMassMatrix(coef, Me_init); 
        fem.FormSystemMatrix(Me_init, boundary_dofs, Mmate); 
        
        Me_solver.iterative_mode = false;
        Me_prec.SetType(mfem::HypreSmoother::Jacobi); 
        fem.SolverConditions(Me_solver, Me_prec); 

        fem.InitializeStiffnessMatrix(cDe, Ke2); 

        PeR = psx;
        PeR.Neg();
        m_nbcCoef = std::make_unique<mfem::ProductCoefficient>(matCoef_R, nbcCoef);

        Be_init = std::make_unique<mfem::ParLinearForm>(fespace.get());
        Be_init->AddDomainIntegrator(new mfem::DomainLFIntegrator(cAe));
        Be_init->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(*m_nbcCoef), nbc_bdr);
        Be_init->Assemble();
        Fet = *Be_init;

        fem.FormLinearSystem(Ke2, boundary_dofs, Cn, Fet, Kmate, X1v, Feb);
    }

    else if(mode_ == sim::CellMode::FULL){
        if (mfem::Mpi::WorldRank() == 0) {std::cout << "Initializing CnE for Full Cell" << std::endl;}

        utils.SetInitialValue(Cn, initial_value);

        // Assemble mass matrix
        mfem::GridFunctionCoefficient coef(&psx);
        fem.InitializeMassMatrix(coef, Me_init); 
        fem.FormSystemMatrix(Me_init, boundary_dofs, Mmate); 
        
        // Configure solver
        Me_solver.iterative_mode = false; 
        Me_prec.SetType(mfem::HypreSmoother::Jacobi); 
        fem.SolverConditions(Me_solver, Me_prec); 

        // Initialize stiffness operator (diffusivity)
        fem.InitializeStiffnessMatrix(cDe, Ke2); 

        // Build force term (domain + boundary contributions)
        Be_init = std::make_unique<mfem::ParLinearForm>(fespace.get());
        Be_init->AddDomainIntegrator(new mfem::DomainLFIntegrator(cAe));
        Be_init->Assemble(); 
        Fet = *Be_init;

        fem.FormLinearSystem(Ke2, boundary_dofs, Cn, Fet, Kmate, X1v, Feb); 
    }
}

void ElectrolyteDiffusion::UpdateConcentration(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx,
                        double gtPsx, mfem::ParGridFunction &weight_elec, const std::vector<ConcentrationBase::PairCoupling> &pair_terms)
{
    (void)weight_elec;
    (void)pair_terms;

    utils.InitializeReaction(Rx, Rxe, (-1.0 * Constants::t_minus));
    cAe.SetGridFunction(&Rxe);
    utils.CalculateReactionInfx(Rxe, infx);

    nbcCoef.constant = infx;
    mfem::ProductCoefficient m_nbcCoef(matCoef_R, nbcCoef);

    Be_init = std::make_unique<mfem::ParLinearForm>(fespace.get());
    Be_init->AddDomainIntegrator(new mfem::DomainLFIntegrator(cAe));
    Be_init->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(m_nbcCoef), nbc_bdr);

    Be_init->Assemble();
    Fet = *Be_init; 


    for (int vi = 0; vi < nV; vi++){
        const double cn_val = Cn(vi);
        De(vi) = psx(vi) * MaterialProperties::Diffusivity(material, cn_val);
    }
    cDe.SetGridFunction(&De);

    fem.Update(Ke2); 
    fem.FormLinearSystem(Ke2, boundary_dofs, Cn, Fet, Kmate, X1v, Feb); 
    Feb *= cfg.dt;

    TmatR.reset(Add(1.0, Mmate, -0.5 * cfg.dt, Kmate));
    TmatL.reset(Add(1.0, Mmate,  0.5 * cfg.dt, Kmate));

    Cn.GetTrueDofs(CeV0);
    TmatR->Mult(CeV0, RHCe);
    RHCe += Feb; 

    Me_solver.SetOperator(*TmatL);
    Me_solver.SetPreconditioner(Me_prec);
    Me_solver.Mult(RHCe, CeVn) ;

    for (int p = 0; p < CeVn.Size(); p++){
        if (CeVn(p) < 0.0){
            (CeVn)(p) = 1.0e-6;} 
    }

    Cn.Distribute(CeVn);
}

void ElectrolyteDiffusion::SaltConservation(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx) {

    CeC = 0.0;
    
    mfem::ParGridFunction CeT(fespace.get());
    
    CeT = Cn;
    CeT *= psx;

    for (int ei = 0; ei < nE; ei++){
        CeT.GetNodalValues(ei,VtxVal); 
        double val = 0.0;
        for (int vt = 0; vt < nC; vt++){val += VtxVal[vt];}        
        EAvg(ei) = val/nC; 
        CeC += EAvg(ei)*EVol(ei); 
    }
    
    MPI_Allreduce(&CeC, &gCeC, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);			
    
    CeAvg = gCeC/gtPse;	

    Cn -= (CeAvg-Ce0);
    MPI_Barrier(MPI_COMM_WORLD);
}