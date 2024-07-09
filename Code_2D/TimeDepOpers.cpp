#include "TimeDepOpers.hpp"
#include "mfem.hpp"

using namespace mfem;
using namespace std;


/*
 * Define a new operator inherited from TimeDependentOperator.  
 * We need this in order to use MFEM's built-in time-dependent ODE solvers.
 * This will be used to solve the concentration evolution in the 
 * cathode and electrolyte, which follow a second order PDE.
 *
 * The general ODE can be written as M du/dt = K u + b
 * An explicit step will solve du/dt = M^{-1}(K u + b)
 * An implicit step will solve du/dt = M^{-1}(K (u+du/dt) + b)
 *                                   = (M - dt K)^{-1}(K u + b)
 * See MFEM examples 9 and 16.
 * 
 * NOTE: We might need a different operator for the evolution in the anode,
 *       which is more complicated and follows a fourth order PDE.
 */

ConductionOperator::ConductionOperator(ParGridFunction &ps, HypreParMatrix &K, HypreParVector &Fb)
   : TimeDependentOperator(K.Height(), K.Width(), (double) 0.0),
     M_solver(K.GetComm()), T_solver(K.GetComm())
{
   const double rel_tol = 1e-8;

   Kmat = K;
   b = Fb;
   z = Fb.CreateCompatibleVector();
   
   GridFunctionCoefficient cp(&ps);
   
   M = new ParBilinearForm(ps.ParFESpace());
   M->AddDomainIntegrator(new MassIntegrator(cp));
   M->Assemble(); // keep sparsity pattern of M and K the same
   M->FormSystemMatrix(ess_tdof_list, Mmat);
   
   M_solver.iterative_mode = false;
   M_solver.SetRelTol(rel_tol);
   M_solver.SetAbsTol(0.0);
   M_solver.SetMaxIter(100);
   M_solver.SetPrintLevel(0);
   M_prec.SetType(HypreSmoother::Jacobi);
   M_solver.SetPreconditioner(M_prec);
   M_solver.SetOperator(Mmat);
   

   T_solver.iterative_mode = false;
   T_solver.SetRelTol(rel_tol);
   T_solver.SetAbsTol(0.0);
   T_solver.SetMaxIter(100);
   T_solver.SetPrintLevel(0);
   T_prec.SetType(HypreSmoother::Jacobi);
   T_solver.SetPreconditioner(T_prec);

}

ConductionOperator::ConductionOperator(ParGridFunction &ps, HypreParMatrix &K, HypreParVector &Fb, Array<int> &ess_list)
   : TimeDependentOperator(K.Height(), K.Width(), (double) 0.0),
     M_solver(K.GetComm()), T_solver(K.GetComm())
{
   const double rel_tol = 1e-8;

   //ess_tdof_list = ess_list;

   Kmat = K;
   b = Fb;
   z = Fb.CreateCompatibleVector();
   
   GridFunctionCoefficient cp(&ps);
   
   M = new ParBilinearForm(ps.ParFESpace());
   M->AddDomainIntegrator(new MassIntegrator(cp));
   M->Assemble(); // keep sparsity pattern of M and K the same
   //M->FormSystemMatrix(ess_tdof_list, Mmat);
   M->FormSystemMatrix(ess_list, Mmat);

   M_solver.iterative_mode = false;
   M_solver.SetRelTol(rel_tol);
   M_solver.SetAbsTol(0.0);
   M_solver.SetMaxIter(100);
   M_solver.SetPrintLevel(0);
   M_prec.SetType(HypreSmoother::Jacobi);
   M_solver.SetPreconditioner(M_prec);
   M_solver.SetOperator(Mmat);
   

   T_solver.iterative_mode = false;
   T_solver.SetRelTol(rel_tol);
   T_solver.SetAbsTol(0.0);
   T_solver.SetMaxIter(100);
   T_solver.SetPrintLevel(0);
   T_prec.SetType(HypreSmoother::Jacobi);
   T_solver.SetPreconditioner(T_prec);

}

void ConductionOperator::Mult(const Vector &u, Vector &du_dt) const
{
   // Compute:
   //    du_dt = M^{-1}*-Ku
   // for du_dt, where K is linearized by using u from the previous timestep
   Kmat.Mult(u, z);
   z.Neg(); // z = -z
   //cout << "Max z: " << z.Max() << endl;
   //cout << "Min z: " << z.Min() << endl;
   z += b;
   //cout << "Max z: " << z.Max() << endl;
   //cout << "Min z: " << z.Min() << endl;
   M_solver.Mult(z, du_dt);
   //cout << "Max dudt: " << du_dt.Max() << endl;
   //cout << "Min duxt: " << du_dt.Min() << endl;
}

void ConductionOperator::ImplicitSolve(const double dt,
                                       const Vector &u, Vector &du_dt)
{
   // Solve the equation:
   //    du_dt = M^{-1}*[-K(u + dt*du_dt)]
   // for du_dt, where K is linearized by using u from the previous timestep
   T = Add(1.0, Mmat, dt, Kmat);
   //current_dt = dt;
   T_solver.SetOperator(*T);
   //MFEM_VERIFY(dt == current_dt, ""); // SDIRK methods use the same dt
   Kmat.Mult(u, z);
   z.Neg();
   z += b;
   T_solver.Mult(z, du_dt);
}
/*
void ConductionOperator::SetParameters(const Vector &u)
{
   
   ParGridFunction u_alpha_gf(&fespace);
   u_alpha_gf.SetFromTrueDofs(u);
   for (int i = 0; i < u_alpha_gf.Size(); i++)
   {
      u_alpha_gf(i) = kappa + alpha*u_alpha_gf(i);
   }

   delete K;
   K = new ParBilinearForm(&fespace);

   GridFunctionCoefficient u_coeff(&u_alpha_gf);

   K->AddDomainIntegrator(new DiffusionIntegrator(u_coeff));
   K->Assemble(0); // keep sparsity pattern of M and K the same
   K->FormSystemMatrix(ess_tdof_list, Kmat);
   delete T;
   T = NULL; // re-compute T on the next ImplicitSolve
   
}
*/

void ConductionOperator::UpdateParams(const HypreParMatrix &K, const HypreParVector &Fb)
{
   Kmat = K;
   b = Fb;
}

ConductionOperator::~ConductionOperator()
{
   //cout << "HERE A DELETE" << endl;
   delete T;
   //cout << "HERE B DELETE" << endl;
   delete M;
   //cout << "HERE C DELETE" << endl;
   //delete K;
   //cout << "HERE D DELETE" << endl;
}



// Implementation of class FE_Evolution
FE_Evolution::FE_Evolution(ParBilinearForm &M_, ParBilinearForm &K_,
                           const Vector &b_)
   : TimeDependentOperator(M_.ParFESpace()->GetTrueVSize()), b(b_),
     M_solver(M_.ParFESpace()->GetComm()),
     z(height)
{
   if (M_.GetAssemblyLevel()==AssemblyLevel::LEGACY)
   {
      M.Reset(M_.ParallelAssemble(), true);
      K.Reset(K_.ParallelAssemble(), true);
   }
   else
   {
      M.Reset(&M_, false);
      K.Reset(&K_, false);
   }

   M_solver.SetOperator(*M);

   Array<int> ess_tdof_list;
   if (M_.GetAssemblyLevel()==AssemblyLevel::LEGACY)
   {
      HypreParMatrix &M_mat = *M.As<HypreParMatrix>();
      HypreParMatrix &K_mat = *K.As<HypreParMatrix>();
      HypreSmoother *hypre_prec = new HypreSmoother(M_mat, HypreSmoother::Jacobi);
      M_prec = hypre_prec;

      dg_solver = new DG_Solver(M_mat, K_mat, *M_.FESpace());
   }
   else
   {
      M_prec = new OperatorJacobiSmoother(M_, ess_tdof_list);
      dg_solver = NULL;
   }

   M_solver.SetPreconditioner(*M_prec);
   M_solver.iterative_mode = false;
   M_solver.SetRelTol(1e-9);
   M_solver.SetAbsTol(0.0);
   M_solver.SetMaxIter(100);
   M_solver.SetPrintLevel(0);
}

// Solve the equation:
//    u_t = M^{-1}(Ku + b),
// by solving associated linear system
//    (M - dt*K) d = K*u + b
void FE_Evolution::ImplicitSolve(const real_t dt, const Vector &x, Vector &k)
{
   K->Mult(x, z);
   z += b;
   dg_solver->SetTimeStep(dt);
   dg_solver->Mult(z, k);
}

void FE_Evolution::Mult(const Vector &x, Vector &y) const
{
   // y = M^{-1} (K x + b)
   K->Mult(x, z);
   z += b;
   M_solver.Mult(z, y);
}

FE_Evolution::~FE_Evolution()
{
   delete M_prec;
   delete dg_solver;
}


