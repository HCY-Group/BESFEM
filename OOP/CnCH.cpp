/**
 * @file CnCH.cpp
 * @brief Implementation of the Cahn Hilliard concentration class for battery simulations.
 */

 #include "CnCH.hpp"
 #include "../code/Constants.hpp"
 #include "mfem.hpp"
 #include <optional>
 
 
 using namespace std;
 
 CnCH::CnCH(Initialize_Geometry &geo, Domain_Parameters &para)
     : Concentrations(geo, para), geometry(geo), domain_parameters(para), fespace(geo.parfespace), Mub(fespace.get()), Mob(fespace.get())
     
     {

        // Ensure AvP is not null
        MFEM_VERIFY(para.AvP != nullptr, "AvP pointer in Domain_Parameters is null.");
        AvP = *para.AvP;

        phV0 = mfem::HypreParVector(fespace.get());
        Lp1 = mfem::HypreParVector(fespace.get());
        RHS = mfem::HypreParVector(fespace.get());
        Lp2 = mfem::HypreParVector(fespace.get());
        MuV = mfem::HypreParVector(fespace.get());

        CpV0 = std::shared_ptr<mfem::HypreParVector>(new mfem::HypreParVector(fespace.get()));
        RHCp = std::shared_ptr<mfem::HypreParVector>(new mfem::HypreParVector(fespace.get()));
        CpVn = std::shared_ptr<mfem::HypreParVector>(new mfem::HypreParVector(fespace.get()));

        PsVc = mfem::HypreParVector(fespace.get());

    
        M_phi = std::make_shared<mfem::HypreParMatrix>();
        MCH_solver = std::make_shared<mfem::CGSolver>(MPI_COMM_WORLD);
        MCH_prec.SetType(mfem::HypreSmoother::Jacobi);

        K_phi = std::make_shared<mfem::HypreParMatrix>();
        K_mu = std::make_shared<mfem::HypreParMatrix>();

        std::ifstream myXfile("../C_Li_X_101.txt");
        std::ifstream mydFfile("../C_Li_M6_101.txt");
        std::ifstream myMBfile("../C_Li_Mb5_101.txt");
        
        for (int i = 0; i < 101; i++) myXfile >> X_101(i);
        for (int i = 0; i < 101; i++) mydFfile >> dF_101(i);

        for (int i = 0; i < 101; ++i) {
            myMBfile >> Mb5_101(i);
            Mb5_101(i) *= 1.0e2 * 2.0 / 3.0;  // Fortran scaling
        }

        // for (int i = 0; i < 101; ++i) {
        //     myMBfile >> Mb5_101(i);
        //     Mb5_101(i) *= 1.0e2 * 2.0 / 3.0; // scaling from Fortran
        //     Mob(i) = Mb5_101(i);        
        // }

        // for (int i = 0; i < 101; i++) myMBfile >> Mb5_101(i);

        // // Scale mobility values to match the FORTRAN behavior
        // // DTbl(2,1:101) = DTbl(2,1:101)*1.0d2*2/3	
        // for (int i = 0; i < 101; i++) {
        //     Mb5_101(i) *= 1.0e2 * 2.0 / 3.0;
        // }
 
     }

double CnCH::Calculate_dF(double ph, const mfem::Vector &X_101, const mfem::Vector &dF_101)
    {
        int idx = std::floor(ph / 0.01);
        if (idx < 0) idx = 0;
        if (idx > 99) idx = 99;
        return dF_101(idx) + (ph - X_101(idx)) / 0.01 * (dF_101(idx + 1) - dF_101(idx));
    }

double CnCH::Calculate_Mobility(double cn, const mfem::Vector &X_101, const mfem::Vector &mob_101)
    {
        // int idx = std::floor(cn / 0.01) + 0;  // +0 to match [0,100] index logic
        // if (idx < 0) idx = 0;
        // if (idx > 99) idx = 99;

        // Tbl_Chk_Pln_3DA
        // Fortran-style clamping
        if (cn < 0.0) {
            cn = 1.0e-6;
        } else if (cn > 1.0) {
            cn = 1.0;
        }

        int idx = std::floor(cn / 0.01);
        if (idx < 0) idx = 0;
        if (idx > 99) idx = 99;

        return mob_101(idx) + (cn - X_101(idx)) / 0.01 * (mob_101(idx + 1) - mob_101(idx));
    }
 
 void CnCH::Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx)
    {
        Concentrations::SetInitialConcentration(Cn, initial_value); // initial value: 2.02d-2
        Concentrations::SetUpSolver(psx, M_phi, *MCH_solver, MCH_prec); // sets up Mass Matrix & Conditions 
    
        psx.GetTrueDofs(PsVc); // Analog to PsVc in CnP
        Cn.GetTrueDofs(*CpV0); // Save initial Cn into CpV0 (like CpV0 in CnP)
    }
 
 void CnCH::TimeStep(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx)
 {
    
    mfem::Array<int> boundary_dofs;  // Empty array = natural Neumann BCs


    // Extract current Cn values
    Cn.GetTrueDofs(*CpV0);

    // Cn.GetTrueDofs(phV0); // Extract true degrees of freedom in the potential field 
    mfem::ConstantCoefficient kap(Constants::eps * Constants::eps); // Diffusion coefficient for the Laplacian term
    
    // boundary_dofs.DeleteAll(); // reset to natural boundary

    // stiffness matrix of phi
    SolverSteps::StiffnessMatrix(kap, boundary_dofs, K_phi);

    // K_phi->Mult(phV0, Lp1); // Lp1 = K_phi * phV0, KU=F, -Lap(u) = f
    K_phi->Mult(*CpV0, Lp1);


    for (int i = 0; i < Cn.Size(); i++) {
        double val = Cn(i);  // get value at DOF i
        Mub(i) = Calculate_dF(val, X_101, dF_101);
        // Mub(i) *= psx(i);
    }

    Mub.Save("Results/Mub"); // save mobility values for debugging

    Mub.GetTrueDofs(MuV);
    MuV += Lp1; // mu = mu_b - Lp1; Eq 4 changed to - from +

    for (int i = 0; i < Cn.Size(); i++) {
        double cn_val = Cn(i);
        Mob(i) = Calculate_Mobility(cn_val, X_101, Mb5_101);
        // Mob(i) *= psx(i); // multiply by psi 	pmN(0:ny+1,0:nx+1,0:nz+1) = psP(0:ny+1,0:nx+1,0:nz+1)*Mb(0:ny+1,0:nx+1,0:nz+1)
    }

    Mob.Save("Results/Mob"); // save mobility values for debugging

    mfem::GridFunctionCoefficient Mob_GFC(&Mob); // updated value

    // boundary_dofs.DeleteAll(); // reset to natural boundary

    // stiffness matrix for mu
    SolverSteps::StiffnessMatrix(Mob_GFC, boundary_dofs, K_mu);

    // Laplace of mu
    K_mu->Mult(MuV, Lp2);
    Lp2.Neg(); 

    // we need a HypreParVector for the reaction term. ==> Rxn
    // the forward and backward rate constants come from tabulating.
	Lp2 += Rx;
	Lp2 *= Constants::dt;

	// right hand side; Eq 5
	// M_phi->Mult(phV0, RHS);
    M_phi->Mult(*CpV0, *RHCp);

    for (int i = 0; i < CpV0->Size(); i++) {
        if (PsVc(i) > 1.0e-5) {
            double rhs_term = Lp2(i) + Rx(i) / Constants::rho;
            rhs_term *= Constants::dt;
            (*RHCp)(i) += rhs_term / PsVc(i);  // SBM normalization
        }
    }

    // (*RHCp) += Lp2;

    // RHS += Lp2; // RHS = M_phi * phV0 + Lp2
  
    // time update
    // MCH_solver->Mult(RHS, phV0); // solve M_phi·phV0^{n+1} = RHS, A*x=B, phV0 = M_phi^-1 * RHS
    MCH_solver->Mult(*RHCp, *CpVn);

    for (int i = 0; i < CpVn->Size(); i++) {
        if (PsVc(i) < 1.0e-5) {
            (*CpVn)(i) = 0.02;
        }
    }

    Cn.Distribute(CpVn.get());


    // Cn.Distribute(phV0);

    // // try to move this outside of the loop -- into main after simulation
    // for (int i = 0; i < Cn.Size(); i++) {
    //     if (psx(i) < 1.0e-5) {
    //         Cn(i) = 0.02;
    //     }
    // }

    // Degree of Lithiation
    Concentrations::LithiationCalculation(Cn, psx);
 
 }
 