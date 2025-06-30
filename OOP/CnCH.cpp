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
     : Concentrations(geo, para), geometry(geo), domain_parameters(para), fespace(geo.parfespace), Mub(fespace.get()), Mob(fespace.get()), Rxc(fespace.get())
     
     {

        Lp1 = mfem::HypreParVector(fespace.get());
        Lp2 = mfem::HypreParVector(fespace.get());
        MuV = mfem::HypreParVector(fespace.get());

        CpV0 = std::shared_ptr<mfem::HypreParVector>(new mfem::HypreParVector(fespace.get()));
        RHCp = std::shared_ptr<mfem::HypreParVector>(new mfem::HypreParVector(fespace.get()));

        PsVc = mfem::HypreParVector(fespace.get());
    
        M_init = std::make_unique<mfem::ParBilinearForm>(fespace.get());
        MmatCH = std::make_shared<mfem::HypreParMatrix>();

        MCH_solver = std::make_shared<mfem::CGSolver>(MPI_COMM_WORLD);
        MCH_prec.SetType(mfem::HypreSmoother::Jacobi);
        
        Grad_MForm = std::make_unique<mfem::ParBilinearForm>(fespace.get());
        Grad_MM = std::make_shared<mfem::HypreParMatrix>();

        B_init = std::make_unique<mfem::ParLinearForm>(fespace.get());

        Grad_EForm = std::make_unique<mfem::ParBilinearForm>(fespace.get());
        Grad_EM = std::make_shared<mfem::HypreParMatrix>();

        Mob = mfem::ParGridFunction(fespace.get());
        Mub = mfem::ParGridFunction(fespace.get());
        Rxc = mfem::ParGridFunction(fespace.get());
        Fct = mfem::ParLinearForm(fespace.get());
        Fcb = mfem::HypreParVector(fespace.get());

        cDp = mfem::GridFunctionCoefficient(&Mob);
        cAp = mfem::GridFunctionCoefficient(&Rxc);

        std::ifstream myXfile("../Inputs/C_Li_X_101.txt"); // ticks
        std::ifstream mydFfile("../Inputs/C_Li_M6_101.txt"); // chemical potential
        std::ifstream myMBfile("../Inputs/C_Li_Mb5_101.txt"); // mobility
        std::ifstream myOCVfile("../Inputs/C_Li_O3_101.txt"); // OCV
        std::ifstream myi0file("../Inputs/C_Li_J2_101.txt"); // i0
        
        for (int i = 0; i < 101; i++) myXfile >> Ticks(i);
        for (int i = 0; i < 101; i++) mydFfile >> chmPot(i);
        for (int i = 0; i < 101; i++) myMBfile >> Mobility(i);
        for (int i = 0; i < 101; i++) myOCVfile >> OCV(i);
        for (int i = 0; i < 101; i++) myi0file >> i0(i);
 
     }

double CnCH::GetTableValues(double cn, const mfem::Vector &ticks, const mfem::Vector &data)
    {
        if (cn < 1.0e-6) {
            cn = 1.0e-6;
        } else if (cn > 1.0) {
            cn = 1.0;
        }

        int idx = std::floor(cn / 0.01);
        if (idx < 0) idx = 0;
        if (idx > 99) idx = 99;

        return data(idx) + (cn - ticks(idx)) / 0.01 * (data(idx + 1) - data(idx));
    }
 
 void CnCH::Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx)
    {
        Concentrations::SetInitialConcentration(Cn, initial_value); // initial value: 2.02d-2
        Concentrations::LithiationCalculation(Cn, psx); // Calculate the initial concentration based on the potential field

        mfem::GridFunctionCoefficient coef(&psx); 
        SolverSteps::InitializeMassMatrix(coef, M_init); 
        SolverSteps::FormSystemMatrix(M_init, boundary_dofs, MmatCH); // Form the system matrix from the bilinear form
        SolverSteps::SolverConditions(MmatCH, *MCH_solver, MCH_prec); // Set up the solver conditions for the mass matrix

        Mob = 1.0e-12; // Initialize mobility to a small value
        Mob *= psx; // Scale

        SolverSteps::InitializeStiffnessMatrix(cDp, Grad_MForm); // Initialize stiffness form for mobility

        mfem::ConstantCoefficient varE(6.7600e-10);
        SolverSteps::InitializeStiffnessMatrix(varE, Grad_EForm); // Initialize stiffness form for energy

        SolverSteps::InitializeForceTerm(cAp, B_init);
        Fct = std::move(*B_init);

        SolverSteps::FormLinearSystem(Grad_EForm, boundary_dofs, Cn, Fct, Grad_EM, X1v, Fcb); // Form the linear system for energy
        SolverSteps::FormLinearSystem(Grad_MForm, boundary_dofs, Mub, Fct, Grad_MM, X1v, Fcb); // Form the linear system for mobility

        psx.GetTrueDofs(PsVc); 
    }
 
 void CnCH::TimeStep(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx)
    {
        Concentrations::CreateReaction(Rx, Rxc, (1.0/Constants::rho));
        cAp.SetGridFunction(&Rxc); // Set the reaction term coefficient for the force term

        SolverSteps::Update(B_init); // Update the force term with the current reaction term
        Fct = std::move(*B_init); // Move the updated force term to Fct

        // Tabulate the chemical potential and mobility values
        for (int i = 0; i < Cn.Size(); i++) {
            double val = Cn(i);  // get value at DOF i
            Mub(i) = GetTableValues(val, Ticks, chmPot);
        }

        for (int i = 0; i < Cn.Size(); i++) {
            double cn_val = Cn(i);
            Mob(i) = psx(i) * GetTableValues(cn_val, Ticks, Mobility);
        }

        Cn.GetTrueDofs(*CpV0);

        // Lap phi
        Grad_EM->Mult(*CpV0, Lp1); // Lp1 = Grad_EM * CpV0
        Mub.GetTrueDofs(MuV); // Get the true degrees of freedom
        MuV += Lp1; // Update MuV with the Laplacian of phi

        SolverSteps::Update(Grad_MForm);
        SolverSteps::FormLinearSystem(Grad_MForm, boundary_dofs, Mub, Fct, Grad_MM, X1v, Fcb); // Form the linear system for updated chemical potential

        Grad_MM->Mult(MuV, Lp2); // Lp2 = Grad_MM * MuV
        Lp2.Neg(); // Negate Lp2 to account for the negative Lap
        Lp2 *= Constants::dt; // Scale Lp2 by the time step

        // Add the reaction term to the right-hand side vector
        Fcb *= Constants::dt;
        Lp2 += Fcb; // Add the free energy contributions to Lp2

        // Update the right-hand side vector for the Cahn-Hilliard equation
        MmatCH->Mult(*CpV0, *RHCp); // Multiply the mass matrix with the current concentration values
        *RHCp += Lp2; // Add the reaction term to the right-hand side vector

        MCH_solver->Mult(*RHCp, *CpV0); // Solve the system M·CpV0 = RHCp

        // Ensure that the concentration values are within the valid range
        for (int i = 0; i < CpV0->Size(); i++) {
            if (PsVc(i) < 1.0e-5) {
                (*CpV0)(i) = 0.0202;
            }
        }

        // update only the solid region concentrations
        for (int i = 0; i < CpV0->Size(); ++i) {
            if ((*CpV0)(i) < 0.0) {
                (*CpV0)(i) = 0.0;
            } else if ((*CpV0)(i) > 1.0) {
                (*CpV0)(i) = 1.0;
            }
        }

        // recover the GridFunction from the HypreParVector
        Cn.Distribute(CpV0.get()); 

        Concentrations::LithiationCalculation(Cn, psx); // Update the degree of lithiation based on the new concentration values
    }
