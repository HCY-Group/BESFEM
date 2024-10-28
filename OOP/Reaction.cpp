#include "Mesh_Handler.hpp"
#include "Constants.hpp"
#include "Concentrations.hpp"
#include "Reaction.hpp"
#include <fstream>
#include <iostream>

using namespace mfem;
using namespace std;

Reaction::Reaction(mfem::ParFiniteElementSpace *fe, MeshHandler &mh, Concentrations &con) 
    : fespace(fe), mesh_handler(mh), AvP(*mh.GetAvP()), AvB(*mh.GetAvB()), concentrations(con)
    


{
    nV = mesh_handler.GetNV();
    
    Rxn = new mfem::ParGridFunction(fespace);

    Dmp = new mfem::ParGridFunction(fespace); // D_minus_plus
    kpl = new mfem::ParGridFunction(fespace); // electrolyte conductivity

    Kl1 = new mfem::ParBilinearForm(fespace);
    B1t = new mfem::ParLinearForm(fespace);
    X1v = new mfem::HypreParVector(fespace);
    B1v = new mfem::HypreParVector(fespace);

    LpCe = new mfem::HypreParVector(fespace);

    Kl2 = new mfem::ParBilinearForm(fespace);

    kap = new mfem::ParGridFunction(fespace); //conductivity in particle

    Kp2 = new mfem::ParBilinearForm(fespace);

    i0C = new mfem::ParGridFunction(fespace);
    OCV = new mfem::ParGridFunction(fespace);
    Kfw = new mfem::ParGridFunction(fespace);
    Kbw = new mfem::ParGridFunction(fespace);


}


void Reaction::Initialize(){

    AvP_PGF = new ParGridFunction(fespace);	
	*AvP_PGF = AvP;

    AvB_PGF = new ParGridFunction(fespace);	
	*AvB_PGF = AvB;

    CreateRx(*Rxn, 0.0);
    SetAvP(*Rxn, *AvP_PGF, 1.0e-08);

    SetZero(*i0C);
    SetZero(*OCV);
    SetZero(*Kfw);
    SetZero(*Kbw);

}

void Reaction::TimeStep() {

    mfem::ParGridFunction &CnP = *concentrations.CnP; // referencing CnP from concentrations class
    mfem::ParGridFunction &CnE = *concentrations.CnE; // referencing CnE from concentrations class

    // std::cout << CnE << std::endl; // CnE is correct here


    mfem::ParGridFunction &psi = concentrations.psi; // referencing psi from concentrations class
    mfem::ParGridFunction &pse = concentrations.pse; // referencing pse from concentrations class

    ElectrolyteConductivity(CnE, pse);

    // std::cout << "Dmp = " << (*Dmp) << std::endl;

    mfem::GridFunctionCoefficient cDm(Dmp);
    mfem::GridFunctionCoefficient cKe(kpl);

    mfem::ParGridFunction &phP = *concentrations.potentials.phP; // referencing phP from potentials class
    mfem::ParGridFunction &phE = *concentrations.potentials.phE; // referencing phE from potentials class

    KMatrix(*Kl1, cDm, boundary_dofs, phE, *B1t, Kdm, *X1v, *B1v);

    mfem::HypreParVector &CeVn = *concentrations.CeVn;
    CnE.GetTrueDofs(CeVn);
    Kdm.Mult(CeVn, *LpCe);

    mfem::Array<int> ess_tdof_list_w = mesh_handler.ess_tdof_list_w;
    KMatrix(*Kl2, cKe, ess_tdof_list_w, phE, *B1t, KmE, *X1v, *B1v);

    mfem::CGSolver &cgPE_solver = *concentrations.potentials.cgPE_solver;
    PCG_Solver(Mpe, cgPE_solver, KmE);

    ParticleConductivity(CnP, psi);

    // std::cout << "kap = " << (*kap) << std::endl;

    mfem::GridFunctionCoefficient cKp(kap);

    mfem::Array<int> ess_tdof_list_e = mesh_handler.ess_tdof_list_e;
    KMatrix(*Kp2, cKp, ess_tdof_list_e, phP, *B1t, KmP, *X1v, *B1v);

    mfem::CGSolver &cgPP_solver = *concentrations.potentials.cgPP_solver;
    PCG_Solver(Mpp, cgPP_solver, KmP);

    ExchangeCurrentDensity(*AvB_PGF, CnP);

    // std::cout << "Kbw = " << (*Kbw) << std::endl;


}


void Reaction::ElectrolyteConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx) {
    
    for (int vi = 0; vi < nV; vi++){

        dffe = exp(-7.02 - 830 * Cn(vi) + 50000 * Cn(vi) * Cn(vi));
        (*Dmp)(vi) = psx(vi) * tc1 * Constants::D0 * dffe;
        (*kpl)(vi) = psx(vi) * tc2 * Constants::D0 * dffe * Cn(vi);

    }


}

void Reaction::KMatrix(mfem::ParBilinearForm &K, mfem::GridFunctionCoefficient &gfc, mfem::Array<int> boundary, mfem::ParGridFunction &potential, 
mfem::ParLinearForm &plf_B, mfem::HypreParMatrix &matrix, mfem::HypreParVector &hpv_X, mfem::HypreParVector &hpv_B){

    K.Update();
    K.AddDomainIntegrator(new DiffusionIntegrator(gfc));
    K.Assemble();
    K.FormLinearSystem(boundary, potential, plf_B, matrix, hpv_X, hpv_B);

}

void Reaction::PCG_Solver(mfem::HypreSmoother &smoother, mfem::CGSolver &cg, mfem::HypreParMatrix &KMatrix){

    smoother.SetType(HypreSmoother::Jacobi);
    cg.SetPreconditioner(smoother);
    cg.SetOperator(KMatrix);


}

void Reaction::ParticleConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx){

    for (int vi = 0; vi < nV; vi++){

        (*kap)(vi) = psx(vi) * (0.01929 + 0.7045 * tanh(2.399 * Cn(vi)) - 0.7238 * tanh(2.412 * Cn(vi)) - 4.2106e-6);

    }

}

void Reaction::ExchangeCurrentDensity(mfem::ParGridFunction &Av_pgf, mfem::ParGridFunction &Cn){

    for (int vi = 0; vi < nV; vi++){

        if(Av_pgf(vi) * Constants::dh > 0.0){
            val = -0.2 * (Cn(vi) - 0.37) - 1.559 - 0.9376 * tanh(8.961 * Cn(vi) - 3.195); // check on this!
            (*i0C)(vi) = pow(10.0, val) * 1.0e-3;

            (*OCV)(vi) = 1.095 * Cn(vi) * Cn(vi) - 8.324e-7 * exp(14.31 * Cn(vi)) + 4.692 * exp(-0.5389 * Cn(vi));

            (*Kfw)(vi) = (*i0C)(vi) / (Constants::Frd * 0.001) * exp(Constants::alp * Constants::Cst1 * (*OCV)(vi));
            (*Kbw)(vi) = (*i0C)(vi) / (Constants::Frd * Cn(vi)) * exp(-Constants::alp * Constants::Cst1 * (*OCV)(vi));

        }
    }
}

void Reaction::CreateRx(mfem::ParGridFunction &Rx, double initial_value) {

    Rx = initial_value;

}

void Reaction::SetAvP(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Av, double value) {

    Rx = Av;
    Rx *= value;
    
}

void Reaction::SetZero(mfem::ParGridFunction &pgf){
    
    for (int i = 0; i < pgf.Size(); ++i) {
        pgf(i) = 0.0; 
    }
}

void Reaction::SetPotentials(Potentials *pot_) {
    potentials = pot_;  // Store the Potentials object for later use
}
