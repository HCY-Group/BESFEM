 #include "../include/CnA.hpp"
 #include "../inputs/Constants.hpp"
 #include "mfem.hpp"
 #include <optional>
 
 
 using namespace std;
 
 CnA::CnA(Initialize_Geometry &geo, Domain_Parameters &para)
     : Concentrations(geo, para), geometry(geo), domain_parameters(para), fespace(geo.parfespace), pmesh(geo.parallelMesh.get()),
      // fields
      Mub(fespace.get()), Mob(fespace.get()), RxA(fespace.get()), gtPsi(para.gtPsi), gtPsA(para.gtPsA),
      // vectors
      Lp1(fespace.get()), Lp2(fespace.get()), MuV(fespace.get()),
      PsVc(fespace.get()), CpV0(fespace.get()), RHCp(fespace.get()), CpVn(fespace.get()),
      // forms
      Fct(fespace.get()), Fcb(fespace.get()),
      // coefficients wrap long-lived fields
      cDp(&Mob), cAp(&RxA),
      // solver
      MCH_solver(MPI_COMM_WORLD)   // value type solver

     
     {
        std::cout << "gtPsA before: " << gtPsA << std::endl;
        
        if (gtPsA < 1.0e-200){
            gtPsA = gtPsi;
        }

        std::ifstream myXfile("../inputs/C_Li_X_101.txt"); ///< Concentration ticks
        std::ifstream mydFfile("../inputs/C_Li_M6_101.txt"); ///< Chemical potential
        std::ifstream myMBfile("../inputs/C_Li_Mb5_101.txt"); ///< Mobility
        std::ifstream myOCVfile("../inputs/C_Li_O3_101.txt"); ///< OCV
        std::ifstream myi0file("../inputs/C_Li_J2_101.txt"); ///< Exchange current density
        
        for (int i = 0; i < 101; i++) myXfile >> Ticks(i);
        for (int i = 0; i < 101; i++) mydFfile >> chmPot(i);
        for (int i = 0; i < 101; i++) myMBfile >> Mobility(i);
        for (int i = 0; i < 101; i++) Mobility(i) *= 100.0 * 2.0/3.0;
        for (int i = 0; i < 101; i++) myOCVfile >> OCV(i);
        for (int i = 0; i < 101; i++) myi0file >> i0(i);
 
     }

double CnA::GetTableValues(double cn, const mfem::Vector &ticks, const mfem::Vector &data)
    {
        if (cn < 1.0e-6) cn = 1.0e-6;
        if (cn > 0.999999) cn = 0.999999;

        int idx = std::floor(cn / 0.01);
        if (idx < 0) idx = 0;
        if (idx > 99) idx = 99;

        return data(idx) + (cn - ticks(idx)) / 0.01 * (data(idx + 1) - data(idx));
    }
 
 void CnA::Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx)
    {
        Concentrations::SetInitialConcentration(Cn, initial_value); // initial value: 2.02d-2
        Concentrations::LithiationCalculation(Cn, psx, gtPsA); // Calculate the initial concentration based on the potential field

        mfem::GridFunctionCoefficient coef(&psx); 
        SolverSteps::InitializeMassMatrix(coef, M_init); 
        SolverSteps::FormSystemMatrix(M_init, boundary_dofs, MmatCH); // Form the system matrix from the bilinear form
        
        MCH_solver.iterative_mode = false; // Enable iterative mode for the solver
        MCH_prec.SetType(mfem::HypreSmoother::Jacobi); // Configure the preconditioner using a Jacobi smoother
        SolverSteps::SolverConditions(MmatCH, MCH_solver, MCH_prec); // Set up the solver conditions for the mass matrix

        SolverSteps::InitializeStiffnessMatrix(cDp, Grad_MForm); // Initialize stiffness form for mobility

        mfem::ConstantCoefficient varE(Constants::gc/pow(Constants::dh, pmesh->Dimension())); // dx2 is in M matrix in MFEM, not in K matrix
        SolverSteps::InitializeStiffnessMatrix(varE, Grad_EForm); // Initialize stiffness form for energy

        SolverSteps::InitializeForceTerm(cAp, B_init);
        Fct = *B_init;

        SolverSteps::FormLinearSystem(Grad_EForm, boundary_dofs, Cn, Fct, Grad_EM, X1v, Fcb); // Form the linear system for energy
        SolverSteps::FormLinearSystem(Grad_MForm, boundary_dofs, Mub, Fct, Grad_MM, X1v, Fcb); // Form the linear system for mobility

        psx.GetTrueDofs(PsVc); 

    }
 
 void CnA::TimeStep(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx)
    {
        Concentrations::CreateReaction(Rx, RxA, (1.0/Constants::rho_A));
        cAp.SetGridFunction(&RxA); // Set the reaction term coefficient for the force term

        SolverSteps::Update(B_init); // Update the force term with the current reaction term
        Fct = *B_init; // Move the updated force term to Fct

        // Tabulate the chemical potential and mobility values
        for (int i = 0; i < Cn.Size(); i++) {
            double val = Cn(i);  // get value at DOF i
            Mub(i) = GetTableValues(val, Ticks, chmPot);
        }

        for (int i = 0; i < Cn.Size(); i++) {
            double cn_val = Cn(i);
            Mob(i) = psx(i) * GetTableValues(cn_val, Ticks, Mobility);
        }

        Cn.GetTrueDofs(CpV0);

        // Lap phi
        Grad_EM.Mult(CpV0, Lp1); // Lp1 = Grad_EM * CpV0
        Mub.GetTrueDofs(MuV); // Get the true degrees of freedom
        MuV += Lp1; // Update MuV with the Laplacian of phi

        // Update stiffness form with new mobility
        cDp.SetGridFunction(&Mob); // Set the mobility coefficient for the stiffness matrix
        SolverSteps::Update(Grad_MForm);
        SolverSteps::FormLinearSystem(Grad_MForm, boundary_dofs, Mub, Fct, Grad_MM, X1v, Fcb); // Form the linear system for updated chemical potential

        Grad_MM.Mult(MuV, Lp2); // Lp2 = Grad_MM * MuV
        Lp2.Neg(); // Negate Lp2 to account for the negative Lap
        Lp2 *= Constants::dt; // Scale Lp2 by the time step

        // Add the reaction term to the right-hand side vector
        Fcb *= Constants::dt;
        Lp2 += Fcb; // Add the free energy contributions to Lp2

        // Update the right-hand side vector for the Cahn-Hilliard equation
        MmatCH.Mult(CpV0, RHCp); // Multiply the mass matrix with the current concentration values
        RHCp += Lp2; // Add the reaction term to the right-hand side vector

        MCH_solver.Mult(RHCp, CpV0); // Solve the system M·CpV0 = RHCp

        // Ensure that the concentration values are within the valid range
        for (int i = 0; i < CpV0.Size(); i++) {
            if (PsVc(i) < 1.0e-5) {
                (CpV0)(i) = Constants::init_CnA; // Reset to initial concentration if potential is too low
            }
        }

        // recover the GridFunction from the HypreParVector
        Cn.Distribute(CpV0); 

        Concentrations::LithiationCalculation(Cn, psx, gtPsA); // Update the degree of lithiation based on the new concentration values
    }
