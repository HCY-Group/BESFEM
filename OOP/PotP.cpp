#include "PotP.hpp"
#include "mfem.hpp"

double BvP = 0.0;

PotP::PotP(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh)
    : Potentials(pm, fe, mh), dbc_e_bdr(pmesh->bdr_attributes.Max()), gtPsi(mesh_handler.GetTotalPsi())
    
    {

    cgPP_solver = new mfem::CGSolver(MPI_COMM_WORLD);
    RpP = new ParGridFunction(fespace); // reaction rate for particle;

    Fpb = mfem::HypreParVector(fespace);
    kap = new mfem::ParGridFunction(fespace); // particle conductivity
    // Kp2 = new mfem::ParBilinearForm(fespace);
    Kp2 = std::make_unique<mfem::ParBilinearForm>(fespace); // Use make_unique

    B1t = mfem::ParLinearForm(fespace);
    X1v = mfem::HypreParVector(fespace);
    B1v = mfem::HypreParVector(fespace);

    dbc_e_bdr.SetSize(pmesh->bdr_attributes.Max());  // Resize based on the mesh
    dbc_e_bdr = 0; // fix this
    dbc_e_bdr[2] = 1; // Applying Dirichlet BC to the east boundary

    mfem:: Array<int> ess_tdof_list_e(0);


    }

mfem::CGSolver* PotP::cgPP_solver = nullptr; // static variable to be used in reaction


void PotP::Initialize(mfem::ParGridFunction &ph, double initial_value)

{
    
    BvP = initial_value;

    // std::cout << BvP << std::endl;
    
    Potentials::SetInitialPotentials(ph, BvP);
    Potentials::SetUpSolver(*cgPP_solver, 1e-7, 82);

    fespace->GetEssentialTrueDofs(dbc_e_bdr, ess_tdof_list_e); // fix this 


    
}


void PotP::TimeStep(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &phx)


{
    ParticleConductivity(Cn, psx);

    mfem::GridFunctionCoefficient cKp(kap);

    Potentials::ImplementBoundaryConditions(dbc_e_Coef, BvP, phx, dbc_e_bdr);

    Kp2 = std::make_unique<ParBilinearForm>(fespace);  // Initialize as a member variable

    Potentials::KMatrix(*Kp2, cKp, ess_tdof_list_e, phx, B1t, KmP, X1v, B1v);
    Potentials::PCG_Solver(Mpp, *cgPP_solver, KmP);

}



void PotP::ParticleConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx){

    for (int vi = 0; vi < nV; vi++){

        (*kap)(vi) = psx(vi) * (0.01929 + 0.7045 * tanh(2.399 * Cn(vi)) - 0.7238 * tanh(2.412 * Cn(vi)) - 4.2106e-6);

    }

}

void PotP::CalculateGlobalError(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror) 

{
    // std::cout << "Entering CalculateGlobalError for PotP" << std::endl;

    Potentials::CreateReaction(Rx, *RpP, Constants::Frd);
    Potentials::ForceTerm(*RpP, ftPotP);
    Potentials::ForceVector(*Kp2, ess_tdof_list_e, phx, ftPotP, KmP, X1v, Fpb, dbc_e_Coef, dbc_e_bdr);
    
    // double previous_globalerror = globalerror_P;

    
    
    Potentials::ErrorCalculation(phx, *cgPP_solver, Fpb, psx, error_P, gerror, gtPsi);

    // std::cout << "Previous globalerror_P: " << previous_globalerror
    //           << ", Updated globalerror_P: " << globalerror_P << std::endl;


}



    	// // Get the local data of the HypreParVector
		// double *Fpb_data = Fpb.GetData();

		// // Print each value of the vector
		// int size1 = Fpb.Size();
		// std::cout << "Fpb values in CGE:" << std::endl;
		// for (int i = 0; i < size1; i++) {
		// 	std::cout << Fpb_data[i] << " ";
		// }
		// std::cout << std::endl;

