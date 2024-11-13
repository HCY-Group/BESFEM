#include "PotE.hpp"
#include "mfem.hpp"
#include "CnE.hpp"



PotE::PotE(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh)
    : Potentials(pm, fe, mh)
    
    {

    cgPE_solver = new mfem::CGSolver(MPI_COMM_WORLD);

    Dmp = new mfem::ParGridFunction(fespace); // D_minus_plus
    kpl = new mfem::ParGridFunction(fespace); // electrolyte conductivity

    Kl1 = new mfem::ParBilinearForm(fespace);
    Kl2 = new mfem::ParBilinearForm(fespace);

    B1t = new mfem::ParLinearForm(fespace);
    X1v = new mfem::HypreParVector(fespace);
    B1v = new mfem::HypreParVector(fespace);

    LpCe = new mfem::HypreParVector(fespace);
    CeVn = new mfem::HypreParVector(fespace);


    }

mfem::CGSolver* PotE::cgPE_solver = nullptr; // static variable to be used in reaction


void PotE::Initialize(mfem::ParGridFunction &ph, double initial_value)

{
    Potentials::SetInitialPotentials(ph, initial_value);
    Potentials::SetUpSolver(*cgPE_solver, 1e-7, 80);

    Vcell = BvP - BvE;

    // std::cout << "Vcell: " << Vcell << std::endl;
    
}

void PotE::TimeStep(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &phx)


{
    ElectrolyteConductivity(Cn, psx);

    mfem::GridFunctionCoefficient cDm(Dmp); 
    mfem::GridFunctionCoefficient cKe(kpl); 

    Potentials::KMatrix(*Kl1, cDm, boundary_dofs, phx, *B1t, Kdm, *X1v, *B1v); 

    Cn.GetTrueDofs(*CeVn); 
    Kdm.Mult(*CeVn, *LpCe);  

    mfem::Array<int> ess_tdof_list_w = mesh_handler.ess_tdof_list_w; 
    Potentials::KMatrix(*Kl2, cKe, ess_tdof_list_w, phx, *B1t, KmE, *X1v, *B1v); 
    Potentials::PCG_Solver(Mpe, *cgPE_solver, KmE); 

}

void PotE::ElectrolyteConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx) {
    
    for (int vi = 0; vi < nV; vi++){

        dffe = exp(-7.02 - 830 * Cn(vi) + 50000 * Cn(vi) * Cn(vi));
        (*Dmp)(vi) = psx(vi) * tc1 * Constants::D0 * dffe;
        (*kpl)(vi) = psx(vi) * tc2 * Constants::D0 * dffe * Cn(vi);

    }

}