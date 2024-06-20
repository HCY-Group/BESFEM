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

void ConductionOperator::Mult(const Vector &u, Vector &du_dt) const
{
   // Compute:
   //    du_dt = M^{-1}*-Ku
   // for du_dt, where K is linearized by using u from the previous timestep
   Kmat.Mult(u, z);
   z.Neg(); // z = -z
   z += b;
   M_solver.Mult(z, du_dt);
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
   delete T;
   delete M;
   delete K;
}



