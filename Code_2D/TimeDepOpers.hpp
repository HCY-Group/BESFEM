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

   virtual void Mult(const Vector &u, Vector &du_dt) const;
   
   /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
       This is the only requirement for high-order SDIRK implicit integration.*/
   virtual void ImplicitSolve(const double dt, const Vector &u, Vector &k);

   //Update K matrix and b vector
   void UpdateParams(const HypreParMatrix &K, const HypreParVector &Fb);

   virtual ~ConductionOperator();
};


#endif
