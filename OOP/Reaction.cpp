#include "Reaction.hpp"
#include "Mesh_Handler.hpp"
#include "mfem.hpp"
#include "CnE.hpp"
#include "PotE.hpp"
#include "PotP.hpp"


Reaction::Reaction(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh)
    : pmesh(pm), fespace(fe), mesh_handler(mh), AvP(*mh.GetAvP()), AvB(*mh.GetAvB())

{
    
    nE = mesh_handler.GetNE();
    nC = mesh_handler.GetNC();
    nV = mesh_handler.GetNV();


    // Dmp = new mfem::ParGridFunction(fespace); // D_minus_plus
    // kpl = new mfem::ParGridFunction(fespace); // electrolyte conductivity
    // kap = new mfem::ParGridFunction(fespace); // particle conductivity


    // Kl1 = new mfem::ParBilinearForm(fespace);
    // Kl2 = new mfem::ParBilinearForm(fespace);
    // Kp2 = new mfem::ParBilinearForm(fespace);

    // B1t = new mfem::ParLinearForm(fespace);
    // X1v = new mfem::HypreParVector(fespace);
    // B1v = new mfem::HypreParVector(fespace);
    
    // LpCe = new mfem::HypreParVector(fespace);

    i0C = new mfem::ParGridFunction(fespace); // exchange current density
    OCV = new mfem::ParGridFunction(fespace); // open circuit voltage
    Kfw = new mfem::ParGridFunction(fespace); // forward reaction constant
    Kbw = new mfem::ParGridFunction(fespace); // backward rection constant

    dPHE = new mfem::ParGridFunction(fespace); // voltage drop



}

void Reaction::Initialize(mfem::ParGridFunction &Rx, double initial_value) {

    SetInitialReaction(Rx, initial_value);
    Rx = AvP;
    Rx *= 1.0e-8;

    // Rx.Print();


}

void Reaction::SetInitialReaction(mfem::ParGridFunction &Rx, double initial_value) {
    
    for (int i = 0; i < Rx.Size(); ++i) {
        Rx(i) = initial_value;
    }

}

// void Reaction::TimeStep(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn1, mfem::ParGridFunction &Cn2, mfem::ParGridFunction &psx1, mfem::ParGridFunction &psx2, mfem::ParGridFunction &phx1, mfem::ParGridFunction &phx2) {

//     // ELECTROLYTE CHUNK OF REACTION
    
//     ElectrolyteConductivity(Cn2, psx2); // move to PotE

//     mfem::GridFunctionCoefficient cDm(Dmp); // move to PotE
//     mfem::GridFunctionCoefficient cKe(kpl); // move to PotE

//     KMatrix(*Kl1, cDm, boundary_dofs, phx2, *B1t, Kdm, *X1v, *B1v); // move to PotE

//     mfem::HypreParVector* CeVn = CnE::GetCeVn(); // move to PotE
//     Cn2.GetTrueDofs(*CeVn); // move to PotE
//     Kdm.Mult(*CeVn, *LpCe); // move to PotE

//     mfem::Array<int> ess_tdof_list_w = mesh_handler.ess_tdof_list_w; // move to PotE
//     KMatrix(*Kl2, cKe, ess_tdof_list_w, phx2, *B1t, KmE, *X1v, *B1v); // move to PotE

//     mfem::CGSolver* cgPE_solver = PotE::GetcgPEsolver(); // move to PotE
//     PCG_Solver(Mpe, *cgPE_solver, KmE); // move to PotE

//     // PARTICLE CHUNK OF REACTION

//     ParticleConductivity(Cn1, psx1); // move to PotP

//     mfem::GridFunctionCoefficient cKp(kap); // move to PotP

//     mfem::Array<int> ess_tdof_list_e = mesh_handler.ess_tdof_list_e; // move to PotP
//     KMatrix(*Kp2, cKp, ess_tdof_list_e, phx1, *B1t, KmP, *X1v, *B1v); // move to PotP

//     mfem::CGSolver* cgPP_solver = PotP::GetcgPPsolver(); // move to PotP
//     PCG_Solver(Mpp, *cgPP_solver, KmP); // move to PotP

//     ExchangeCurrentDensity(AvB, Cn1); // Keep in Reaction as its own function

// }

void Reaction::ButlerVolmer(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn1, mfem::ParGridFunction &Cn2, mfem::ParGridFunction &phx1, mfem::ParGridFunction &phx2){

    for (int vi = 0; vi < nV; vi++){
        if ( AvB(vi) * Constants::dh > 0.0 ){

            (*dPHE)(vi) = phx1(vi) - phx2(vi);
            Rx(vi) = AvP(vi) * ((*Kfw)(vi)*Cn2(vi)*exp(-Constants::alp*Constants::Cst1*(*dPHE)(vi)) - \
					                   (*Kbw)(vi)*Cn1(vi)*exp( Constants::alp*Constants::Cst1*(*dPHE)(vi)));

        }
    }

}

// void Reaction::ElectrolyteConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx) {
    
//     for (int vi = 0; vi < nV; vi++){

//         dffe = exp(-7.02 - 830 * Cn(vi) + 50000 * Cn(vi) * Cn(vi));
//         (*Dmp)(vi) = psx(vi) * tc1 * Constants::D0 * dffe;
//         (*kpl)(vi) = psx(vi) * tc2 * Constants::D0 * dffe * Cn(vi);

//     }

// }

// void Reaction::ParticleConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx){

//     for (int vi = 0; vi < nV; vi++){

//         (*kap)(vi) = psx(vi) * (0.01929 + 0.7045 * tanh(2.399 * Cn(vi)) - 0.7238 * tanh(2.412 * Cn(vi)) - 4.2106e-6);

//     }

// }

// void Reaction::KMatrix(mfem::ParBilinearForm &K, mfem::GridFunctionCoefficient &gfc, mfem::Array<int> boundary, mfem::ParGridFunction &potential, 
// mfem::ParLinearForm &plf_B, mfem::HypreParMatrix &matrix, mfem::HypreParVector &hpv_X, mfem::HypreParVector &hpv_B){

//     K.Update();
//     K.AddDomainIntegrator(new DiffusionIntegrator(gfc));
//     K.Assemble();
//     K.FormLinearSystem(boundary, potential, plf_B, matrix, hpv_X, hpv_B);

// }

// void Reaction::PCG_Solver(mfem::HypreSmoother &smoother, mfem::CGSolver &cg, mfem::HypreParMatrix &KMatrix){

//     smoother.SetType(HypreSmoother::Jacobi);
//     cg.SetPreconditioner(smoother);
//     cg.SetOperator(KMatrix);


// }

// rate constants and exchange current density at interface
void Reaction::ExchangeCurrentDensity(mfem::ParGridFunction &Cn){

    for (int vi = 0; vi < nV; vi++){

        if(AvB(vi) * Constants::dh > 0.0){
            double val = -0.2 * (Cn(vi) - 0.37) - 1.559 - 0.9376 * tanh(8.961 * Cn(vi) - 3.195); // check on this!
            (*i0C)(vi) = pow(10.0, val) * 1.0e-3;

            (*OCV)(vi) = 1.095 * Cn(vi) * Cn(vi) - 8.324e-7 * exp(14.31 * Cn(vi)) + 4.692 * exp(-0.5389 * Cn(vi));

            (*Kfw)(vi) = (*i0C)(vi) / (Constants::Frd * 0.001) * exp(Constants::alp * Constants::Cst1 * (*OCV)(vi));
            (*Kbw)(vi) = (*i0C)(vi) / (Constants::Frd * Cn(vi)) * exp(-Constants::alp * Constants::Cst1 * (*OCV)(vi));

        }
    }
}