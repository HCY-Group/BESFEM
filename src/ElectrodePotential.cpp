#include "../include/ElectrodePotential.hpp"
#include "../include/Constants.hpp"
#include "../include/MaterialProperties.hpp"
#include "../include/SimTypes.hpp"
#include "mfem.hpp"

ElectrodePotential::ElectrodePotential(Initialize_Geometry &geo, Domain_Parameters &para, BoundaryConditions &bc, sim::Electrode electrode, sim::MaterialType material, const SimulationConfig &cfg)
    : PotentialBase(geo,para), geometry(geo), domain_parameters(para), boundary_conditions(bc), cfg(cfg), utils(geo,para, cfg), 
    fem(geo.parfespace), electrode(electrode), material(material), kap(fespace.get()), RpP(fespace.get()), pP0(fespace.get())

    
    {
        cgPP_solver = mfem::CGSolver(MPI_COMM_WORLD);
        B1t = mfem::ParLinearForm(fespace.get());
        X1v = mfem::HypreParVector(fespace.get());
        B1v = mfem::HypreParVector(fespace.get());
        Fpb = mfem::HypreParVector(fespace.get());
        Xs0 = mfem::HypreParVector(fespace.get()); 
        RpP = mfem::ParGridFunction(fespace.get());

        Bp2 = std::make_unique<mfem::ParLinearForm>(fespace.get());
        Kp2 = std::make_unique<mfem::ParBilinearForm>(fespace.get()); 

        kap = mfem::ParGridFunction(fespace.get()); 
        cKp = mfem::GridFunctionCoefficient(&kap); 
        cRp = mfem::GridFunctionCoefficient(&RpP);
        Fpt = mfem::ParLinearForm(fespace.get());
        pP0 = mfem::ParGridFunction(fespace.get()); 

        if (electrode == sim::Electrode::ANODE)
        {
            gtPs = para.gtPsiA;
            if (gtPs < 1.0e-200) { gtPs = para.gtPsiA; }
        }
        else if (electrode == sim::Electrode::CATHODE)
        {
            gtPs = para.gtPsiC;
            if (gtPs < 1.0e-200) { gtPs = para.gtPsiC; }
        }
        else
        {
            mfem::mfem_error("ElectrodePotential only supports ANODE or CATHODE.");
        }
    }

void ElectrodePotential::SetupField(mfem::ParGridFunction &ph, double initial_value, mfem::ParGridFunction &psx)
{

    Bv = initial_value;
    utils.SetInitialValue(ph, Bv);

    for (int vi = 0; vi < kap.Size(); vi++)
    {
        kap(vi) = psx(vi) * MaterialProperties::Conductivity(material, 0.0);
    }

    fem.InitializeStiffnessMatrix(cKp, Kp2);
    mfem::ConstantCoefficient dbc_Coef(Bv);

    if (electrode == sim::Electrode::ANODE)
    {
        ph.ProjectBdrCoefficient(dbc_Coef, boundary_conditions.dbc_w_bdr);
        fespace->GetEssentialTrueDofs(boundary_conditions.dbc_w_bdr, ess_tdof_list);
    }
    else
    {
        ph.ProjectBdrCoefficient(dbc_Coef, boundary_conditions.dbc_e_bdr);
        fespace->GetEssentialTrueDofs(boundary_conditions.dbc_e_bdr, ess_tdof_list);
    }

    fem.FormLinearSystem(Kp2, ess_tdof_list, ph, B1t, KmP, X1v, B1v);

    Mpp = std::make_unique<mfem::HypreBoomerAMG>(KmP);
    Mpp->SetPrintLevel(0);
    cgPP_solver.SetRelTol(1e-6);
    cgPP_solver.SetMaxIter(102);
    cgPP_solver.SetPreconditioner(*Mpp);
    cgPP_solver.SetOperator(KmP);

    fem.InitializeForceTerm(cRp, Bp2);
    Bp2->Assemble();
    Fpt = *Bp2;

    fem.FormLinearSystem(Kp2, ess_tdof_list, ph, Fpt, KmP, X1v, Fpb);
}

void ElectrodePotential::AssembleSystem(const std::vector<mfem::ParGridFunction*> &Cn_groups, const std::vector<mfem::ParGridFunction*> &psi_groups, const std::vector<sim::MaterialType> &materials, mfem::ParGridFunction &potential)
{
    mfem::ConstantCoefficient dbc_Coef(Bv);
    cgPP_solver.SetPreconditioner(*Mpp);
    cgPP_solver.SetOperator(KmP);

    ParticleConductivityMulti(Cn_groups, psi_groups, materials);
    fem.Update(Kp2);

    if (electrode == sim::Electrode::ANODE)
    {
        potential.ProjectBdrCoefficient(dbc_Coef, boundary_conditions.dbc_w_bdr);
        fespace->GetEssentialTrueDofs(boundary_conditions.dbc_w_bdr, ess_tdof_list);
    }
    else
    {
        potential.ProjectBdrCoefficient(dbc_Coef, boundary_conditions.dbc_e_bdr);
        fespace->GetEssentialTrueDofs(boundary_conditions.dbc_e_bdr, ess_tdof_list);
    }

    fem.FormLinearSystem(Kp2, ess_tdof_list, potential, B1t, KmP, X1v, B1v);

    Mpp = std::make_unique<mfem::HypreBoomerAMG>(KmP);
    Mpp->SetPrintLevel(0);
    cgPP_solver.SetPreconditioner(*Mpp);
    cgPP_solver.SetOperator(KmP);
}

void ElectrodePotential::ParticleConductivityMulti(const std::vector<mfem::ParGridFunction*> &Cn_groups, const std::vector<mfem::ParGridFunction*> &psi_groups, const std::vector<sim::MaterialType> &materials)
{
    for (int vi = 0; vi < kap.Size(); vi++)
    {
        double kap_sum = 0.0;

        for (size_t j = 0; j < Cn_groups.size(); j++)
        {
            const double c_val = (*Cn_groups[j])(vi);
            const double psi_val = (*psi_groups[j])(vi);

            kap_sum += psi_val * MaterialProperties::Conductivity(materials[j], c_val);
        }

        kap(vi) = kap_sum;
    }
}

void ElectrodePotential::UpdatePotential(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror)
{
    RpP = Rx;
    RpP *= Constants::Frd;

    Bp2->Update();
    Bp2->Assemble();
    Fpt = *Bp2;

    mfem::ConstantCoefficient dbc_Coef(Bv);
    if (electrode == sim::Electrode::ANODE)
    {
        phx.ProjectBdrCoefficient(dbc_Coef, boundary_conditions.dbc_w_bdr);
    }
    else
    {
        phx.ProjectBdrCoefficient(dbc_Coef, boundary_conditions.dbc_e_bdr);
    }

    fem.FormLinearSystem(Kp2, ess_tdof_list, phx, Fpt, KmP, X1v, Fpb);

    pP0 = phx;
    pP0.GetTrueDofs(Xs0);

    cgPP_solver.Mult(Fpb, Xs0);

    phx.Distribute(Xs0);
    utils.CalculateGlobalError(pP0, phx, psx, gerror, gtPs);
}

void ElectrodePotential::AssembleSystem(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &potential)
{
    std::vector<mfem::ParGridFunction*> Cn_groups  = { &Cn };
    std::vector<mfem::ParGridFunction*> psi_groups = { &psx };
    std::vector<sim::MaterialType> materials = { material };

    AssembleSystem(Cn_groups, psi_groups, materials, potential);
}