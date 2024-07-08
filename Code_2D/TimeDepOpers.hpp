#ifndef TIME_DEP_OPERS_HPP
#define TIME_DEP_OPERS_HPP

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include "mpi.h"

#include <cmath>

// #include <stdio.h>

using namespace std;
using namespace mfem;



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
class ConductionOperator : public TimeDependentOperator
{
protected:
   //ParFiniteElementSpace &fespace;
   Array<int> ess_tdof_list; // this list remains empty for pure Neumann b.c.

   ParBilinearForm *M;
   ParBilinearForm *K;

   HypreParVector b;

   HypreParMatrix Mmat;
   HypreParMatrix Kmat;
   HypreParMatrix *T; // T = M + dt K
   //double current_dt;

   CGSolver M_solver;    // Krylov solver for inverting the mass matrix M
   HypreSmoother M_prec; // Preconditioner for the mass matrix M

   CGSolver T_solver;    // Implicit solver for T = M + dt K
   HypreSmoother T_prec; // Preconditioner for the implicit solver

   mutable HypreParVector z; // auxiliary vector

public:
   ConductionOperator(ParGridFunction &ps, HypreParMatrix &K, HypreParVector &Fb);
   ConductionOperator(ParGridFunction &ps, HypreParMatrix &K, HypreParVector &Fb, Array<int> &ess_list);

   virtual void Mult(const Vector &u, Vector &du_dt) const;
   
   /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
       This is the only requirement for high-order SDIRK implicit integration.*/
   virtual void ImplicitSolve(const double dt, const Vector &u, Vector &k);

   //Update K matrix and b vector
   void UpdateParams(const HypreParMatrix &K, const HypreParVector &Fb);

   virtual ~ConductionOperator();
};












// Discontinuous Galerkin Solver used in Advection
class DG_Solver : public Solver
{
private:
   HypreParMatrix &M, &K;
   SparseMatrix M_diag;
   HypreParMatrix *A;
   GMRESSolver linear_solver;
   Solver *prec;
   real_t dt;
public:
   DG_Solver(HypreParMatrix &M_, HypreParMatrix &K_, const FiniteElementSpace &fes)
      : M(M_),
        K(K_),
        A(NULL),
        linear_solver(M.GetComm()),
        dt(-1.0)
   {
      int block_size = fes.GetFE(0)->GetDof();
      //if (prec_type == PrecType::ILU)
      //{
         prec = new BlockILU(block_size,
                             BlockILU::Reordering::MINIMUM_DISCARDED_FILL);
      //}
      //else if (prec_type == PrecType::AIR)
      //{
//#if MFEM_HYPRE_VERSION >= 21800
//         prec = new AIR_prec(block_size);
//#else
//         MFEM_ABORT("Must have MFEM_HYPRE_VERSION >= 21800 to use AIR.\n");
//#endif
//      }
      linear_solver.iterative_mode = false;
      linear_solver.SetRelTol(1e-9);
      linear_solver.SetAbsTol(0.0);
      linear_solver.SetMaxIter(100);
      linear_solver.SetPrintLevel(0);
      linear_solver.SetPreconditioner(*prec);

      M.GetDiag(M_diag);
   }

   void SetTimeStep(real_t dt_)
   {
      if (dt_ != dt)
      {
         dt = dt_;
         // Form operator A = M - dt*K
         delete A;
         A = Add(-dt, K, 0.0, K);
         SparseMatrix A_diag;
         A->GetDiag(A_diag);
         A_diag.Add(1.0, M_diag);
         // this will also call SetOperator on the preconditioner
         linear_solver.SetOperator(*A);
      }
   }

   void SetOperator(const Operator &op)
   {
      linear_solver.SetOperator(op);
   }

   virtual void Mult(const Vector &x, Vector &y) const
   {
      linear_solver.Mult(x, y);
   }

   ~DG_Solver()
   {
      delete prec;
      delete A;
   }
};


// Advection TDO
/** A time-dependent operator for the right-hand side of the ODE. The DG weak
    form of du/dt = -v.grad(u) is M du/dt = K u + b, where M and K are the mass
    and advection matrices, and b describes the flow on the boundary. This can
    be written as a general ODE, du/dt = M^{-1} (K u + b), and this class is
    used to evaluate the right-hand side. */
class FE_Evolution : public TimeDependentOperator
{
private:
   OperatorHandle M, K;
   const Vector &b;
   Solver *M_prec;
   CGSolver M_solver;
   DG_Solver *dg_solver;

   mutable Vector z;

public:
   FE_Evolution(ParBilinearForm &M_, ParBilinearForm &K_, const Vector &b_);

   virtual void Mult(const Vector &x, Vector &y) const;
   virtual void ImplicitSolve(const real_t dt, const Vector &x, Vector &k);

   virtual ~FE_Evolution();
};








#endif

