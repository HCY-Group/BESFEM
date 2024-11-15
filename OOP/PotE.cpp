#include "PotE.hpp"
#include "mfem.hpp"
#include "CnE.hpp"

double BvE = 0.0;

PotE::PotE(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh)
    : Potentials(pm, fe, mh), dbc_w_bdr(pmesh->bdr_attributes.Max())
    
    {

    cgPE_solver = new mfem::CGSolver(MPI_COMM_WORLD);
    RpE = new ParGridFunction(fespace); // reaction rate for particle;

    Flb = mfem::HypreParVector(fespace);

    Dmp = new mfem::ParGridFunction(fespace); // D_minus_plus
    kpl = new mfem::ParGridFunction(fespace); // electrolyte conductivity

    Kl1 = std::make_unique<mfem::ParBilinearForm>(fespace); // Use make_unique
    Kl2 = std::make_unique<mfem::ParBilinearForm>(fespace); // Use make_unique

    B1t = mfem::ParLinearForm(fespace);
    X1v = mfem::HypreParVector(fespace);
    B1v = mfem::HypreParVector(fespace);

    LpCe = new mfem::HypreParVector(fespace);
    CeVn = new mfem::HypreParVector(fespace);

    dbc_w_bdr.SetSize(pmesh->bdr_attributes.Max());  // Resize based on the mesh
    dbc_w_bdr = 0; // fix this
    dbc_w_bdr[0] = 1;

    Array<int> ess_tdof_list_w(0);

    }

mfem::CGSolver* PotE::cgPE_solver = nullptr; // static variable to be used in reaction


void PotE::Initialize(mfem::ParGridFunction &ph, double initial_value)

{
    BvE = initial_value;

    // std::cout << BvE << std::endl;

    
    Potentials::SetInitialPotentials(ph, BvE);
    Potentials::SetUpSolver(*cgPE_solver, 1e-7, 80);

    fespace->GetEssentialTrueDofs(dbc_w_bdr, ess_tdof_list_w); // fix this 

    
}

void PotE::TimeStep(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &phx)


{
    ElectrolyteConductivity(Cn, psx);

    mfem::GridFunctionCoefficient cDm(Dmp); 
    mfem::GridFunctionCoefficient cKe(kpl); 
    
    Kl1 = std::make_unique<ParBilinearForm>(fespace);  // Initialize as a member variable

    Potentials::KMatrix(*Kl1, cDm, boundary_dofs, phx, B1t, Kdm, X1v, B1v); 

    Cn.GetTrueDofs(*CeVn); 
    Kdm.Mult(*CeVn, *LpCe);  

    Kl2 = std::make_unique<ParBilinearForm>(fespace);  // Initialize as a member variable

    Potentials::ImplementBoundaryConditions(dbc_w_Coef, BvE, phx, dbc_w_bdr);
    Potentials::KMatrix(*Kl2, cKe, ess_tdof_list_w, phx, B1t, KmE, X1v, B1v); 
    Potentials::PCG_Solver(Mpe, *cgPE_solver, KmE); 

}

void PotE::ElectrolyteConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx) {
    
    for (int vi = 0; vi < nV; vi++){

        dffe = exp(-7.02 - 830 * Cn(vi) + 50000 * Cn(vi) * Cn(vi));
        (*Dmp)(vi) = psx(vi) * tc1 * Constants::D0 * dffe;
        (*kpl)(vi) = psx(vi) * tc2 * Constants::D0 * dffe * Cn(vi);

    }

}

void PotE::CalculateGlobalError(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx) 

{

    Potentials::CreateReaction(Rx, *RpE, -1.0);
    // RpE->Print(std::cout);

    Potentials::ForceTerm(*RpE, ftPotE);
    Potentials::ForceVector(*Kl2, ess_tdof_list_w, phx, ftPotE, KmE, X1v, Flb, dbc_w_Coef, dbc_w_bdr);

        // 	// Get the local data of the HypreParVector
		// double *Flb_data = Flb.GetData();

		// // Print each value of the vector
		// int size1 = Flb.Size();
		// std::cout << "Flb values in CGE:" << std::endl;
		// for (int i = 0; i < size1; i++) {
		// 	std::cout << Flb_data[i] << " ";
		// }
		// std::cout << std::endl;

}