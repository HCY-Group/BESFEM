#include "Mesh_Handler.hpp"
#include "Constants.hpp"
#include "Concentrations.hpp"
#include "Reaction.hpp"
#include <fstream>
#include <iostream>

using namespace mfem;
using namespace std;

Reaction::Reaction(mfem::ParFiniteElementSpace *fe, MeshHandler &mh, Concentrations &con) 
    : fespace(fe), mesh_handler(mh), AvP(*mh.GetAvP()), concentrations(con)
    


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

}


void Reaction::Initialize(){

    AvP_PGF = new ParGridFunction(fespace);	
	*AvP_PGF = AvP;

    CreateRx(*Rxn, 0.0);
    SetAvP(*Rxn, *AvP_PGF, 1.0e-08);

}

void Reaction::TimeStep() {

    mfem::ParGridFunction &CnE = *concentrations.CnE;
    mfem::ParGridFunction &pse = concentrations.pse;
    ElectrolyteConductivity(CnE, pse);

    mfem::GridFunctionCoefficient cDm(Dmp);
    mfem::GridFunctionCoefficient cKe(kpl);

    mfem::ParGridFunction &phE = *concentrations.potentials.phE;
    KMatrix(*Kl1, cDm, boundary_dofs, phE, *B1t, Kdm, *X1v, *B1v);

    mfem::HypreParVector &CeVn = *concentrations.CeVn;
    CnE.GetTrueDofs(CeVn);
    Kdm.Mult(CeVn, *LpCe); // seg fault




}


void Reaction::ElectrolyteConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx) {

    for (int vi = 0; vi < nV; vi++){

        dffe = exp(-7.02 - 830 * Cn(vi) + 5000 * Cn(vi) * Cn(vi));
        (*Dmp)(vi) = psx(vi) * tc1 * Constants::D0 * dffe;
        (*kpl)(vi) = psx(vi) * tc2 * Constants::D0 * dffe * Cn(vi);
    }

}

void Reaction::KMatrix(mfem::ParBilinearForm &K, mfem::GridFunctionCoefficient &gfc, mfem::Array<int> boundary, mfem::ParGridFunction &potential, 
mfem::ParLinearForm &plf_B, mfem::HypreParMatrix &matrix, mfem::HypreParVector &hpv_X, mfem::HypreParVector &hpv_B){

    HypreParMatrix Khpm;
    Khpm = matrix;

    K.Update();
    K.AddDomainIntegrator(new DiffusionIntegrator(gfc));
    K.Assemble();
    K.FormLinearSystem(boundary, potential, plf_B, Khpm, hpv_X, hpv_B);

}


void Reaction::CreateRx(mfem::ParGridFunction &Rx, double initial_value) {

    Rx = initial_value;

}

void Reaction::SetAvP(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Av, double value) {

    Rx = Av;
    Rx *= value;
    
}

void Reaction::SetPotentials(Potentials *pot_) {
    potentials = pot_;  // Store the Potentials object for later use
}
