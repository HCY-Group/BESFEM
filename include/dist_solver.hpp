// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

/** @file dist_solver.hpp
 * @brief MFEM distance solvers and screened-Poisson filtering utilities.
 */
#ifndef MFEM_DIST_SOLVER_HPP
#define MFEM_DIST_SOLVER_HPP

#include "mfem.hpp"

#ifdef MFEM_USE_MPI

namespace mfem
{

namespace common
{

/** @brief Estimate typical element size from global mean element volume.
 * @param pmesh Parallel mesh; all communicator ranks participate.
 * @return Geometry-dependent length derived from mean element volume.
 */
real_t AvgElementSize(ParMesh &pmesh);

/** @brief Interface for scalar and vector distance-to-level-set calculations. */
class DistanceSolver
{
protected:
   /**
    * @brief Convert a scalar distance field to a vector distance approximation.
    * @param dist_s Scalar distance field.
    * @param[out] dist_v Vector distance field on a compatible space.
    */
   void ScalarDistToVector(ParGridFunction &dist_s, ParGridFunction &dist_v);

public:
   /// 0 = nothing, 1 = main solver only, 2 = full (solver + preconditioner).
   IterativeSolver::PrintLevel print_level; ///< Solver and preconditioner verbosity.

   /**
    * @brief Construct the common distance-solver interface.
    */
   DistanceSolver() { }
   /**
    * @brief Allow cleanup through a distance-solver base pointer.
    */
   virtual ~DistanceSolver() { }

   /// Computes a scalar ParGridFunction which is the length of the shortest path
   /// to the zero level set of the given Coefficient. It is expected that the
   /// given [distance] has a valid (scalar) ParFiniteElementSpace, and that the
   /// result is computed in the same space. Some implementations may output a
   /// "signed" distance, i.e., the distance has different signs on both sides of
   /// the zero level set.
   /**
    * @brief Compute scalar distance to the supplied zero level set.
    * @param zero_level_set Coefficient whose zero contour defines the target.
    * @param[out] distance Result on an initialized scalar parallel space.
    */
   virtual void ComputeScalarDistance(Coefficient &zero_level_set,
                                      ParGridFunction &distance) = 0;

   /// Computes a vector ParGridFunction where the magnitude is the length of the
   /// shortest path to the zero level set of the given Coefficient, and the
   /// direction is the starting direction of the shortest path. It is expected
   /// that the given [distance] has a valid (vector) ParFiniteElementSpace, and
   /// that the result is computed in the same space.
   /**
    * @brief Compute distance magnitude and direction to the zero level set.
    * @param zero_level_set Coefficient whose zero contour defines the target.
    * @param[out] distance Result on an initialized vector parallel space.
    */
   virtual void ComputeVectorDistance(Coefficient &zero_level_set,
                                      ParGridFunction &distance);
};


// K. Crane et al: "Geodesics in Heat: A New Approach to Computing Distance
// Based on Heat Flow", DOI:10.1145/2516971.2516977.
/** @brief Approximate distance using heat diffusion. */
class HeatDistanceSolver : public DistanceSolver
{
public:
   /**
    * @brief Configure a heat-based distance solver.
    * @param diff_coeff Heat diffusion parameter.
    */
   HeatDistanceSolver(real_t diff_coeff)
      : DistanceSolver(), parameter_t(diff_coeff), smooth_steps(0),
        diffuse_iter(1), transform(true), vis_glvis(false) { }

   /// The computed distance is not "signed". In addition to the standard usage
   /// (with zero level sets), this function can be applied to point sources when
   /// transform = false.
   /**
    * @brief Compute unsigned heat-based distance.
    * @param zero_level_set Target level set, or point source when transform is false.
    * @param[out] distance Scalar distance field.
    */
   void ComputeScalarDistance(Coefficient &zero_level_set,
                              ParGridFunction &distance);

   real_t parameter_t; ///< Heat diffusion parameter.
   int smooth_steps; ///< Number of source-smoothing steps.
   int diffuse_iter; ///< Number of heat-diffusion iterations.
   bool transform; ///< Whether to transform the source level set.
   bool vis_glvis; ///< Whether to visualize intermediate fields with GLVis.
};

// A. Belyaev et al: "On Variational and PDE-based Distance Function
// Approximations", Section 6, DOI:10.1111/cgf.12611.
// This solver is computationally cheap, but is accurate for distance
// approximations only near the zero level set.
/** @brief Approximate distance near a zero level set by gradient normalization. */
class NormalizationDistanceSolver : public DistanceSolver
{
private:

   /** @brief Coefficient evaluating the normalized level-set approximation. */
   class NormalizationCoeff : public Coefficient
   {
   private:
      ParGridFunction &u; ///< Borrowed level-set field.

   public:
      /**
       * @brief Borrow the level-set field to normalize.
       * @param u_gf Level-set grid function; must outlive this coefficient.
       */
      NormalizationCoeff(ParGridFunction &u_gf) : u(u_gf) { }
      /**
       * @brief Evaluate the scalar coefficient at an integration point.
       * @param T Element transformation.
       * @param ip Reference integration point.
       * @return Coefficient value.
       */
      real_t Eval(ElementTransformation &T,
                  const IntegrationPoint &ip) override;
   };

public:
   /**
    * @brief Construct the near-interface normalization solver.
    */
   NormalizationDistanceSolver() { }

   /**
    * @brief Compute a normalized near-interface distance approximation.
    * @param u_coeff Input level-set coefficient.
    * @param[out] dist Scalar distance approximation.
    */
   void ComputeScalarDistance(Coefficient& u_coeff, ParGridFunction& dist);
};


// A. Belyaev et al: "On Variational and PDE-based Distance Function
// Approximations", Section 7, DOI:10.1111/cgf.12611.
/** @brief Approximate signed distance using successive p-Laplacian solves. */
class PLapDistanceSolver : public DistanceSolver
{
public:
   /**
    * @brief Set p-Laplacian continuation and Newton solver controls.
    * @param maxp_ Maximum power in the continuation.
    * @param newton_iter_ Maximum Newton iterations per solve.
    * @param rtol Newton relative tolerance.
    * @param atol Newton absolute tolerance.
    */
   PLapDistanceSolver(int maxp_ = 30, int newton_iter_ = 10,
                      real_t rtol = 1e-7, real_t atol = 1e-12)
      : maxp(maxp_), newton_iter(newton_iter_),
        newton_rel_tol(rtol), newton_abs_tol(atol) { }

   /**
    * @brief Set the maximum continuation power.
    * @param new_pp Maximum p value.
    */
   void SetMaxPower(int new_pp) { maxp = new_pp; }

   /// The computed distance is "signed".
   /**
    * @brief Compute signed distance by p-Laplacian continuation.
    * @param func Input level-set coefficient.
    * @param[out] fdist Signed distance approximation.
    */
   void ComputeScalarDistance(Coefficient& func, ParGridFunction& fdist);

private:
   int maxp; ///< Maximum continuation power p.
   const int newton_iter; ///< Newton iteration limit.
   const real_t newton_rel_tol; ///< Newton relative tolerance.
   const real_t newton_abs_tol; ///< Newton absolute tolerance.
};

/** @brief Evaluate the negative normalized gradient of a grid function. */
class NormalizedGradCoefficient : public VectorCoefficient
{
private:
   const GridFunction &u; ///< Borrowed scalar field.

public:
   /**
    * @brief Borrow a field for normalized-gradient evaluation.
    * @param u_gf Scalar field; must outlive this coefficient.
    * @param dim Vector dimension.
    */
   NormalizedGradCoefficient(const GridFunction &u_gf, int dim)
      : VectorCoefficient(dim), u(u_gf) { }

   using VectorCoefficient::Eval;

   /**
    * @brief Evaluate the negative gradient divided by its norm plus 1e-12.
    * @param[out] V Normalized gradient vector.
    * @param T Element transformation.
    * @param ip Reference integration point.
    */
   void Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip)
   {
      T.SetIntPoint(&ip);

      u.GetGradient(T, V);
      const real_t norm = V.Norml2() + 1e-12;
      V /= -norm;
   }
};


// Product of the modulus of the first coefficient and the second coefficient
/** @brief Evaluate the absolute value of one coefficient times another. */
class PProductCoefficient : public Coefficient
{
private:
   Coefficient &basef; ///< Borrowed base coefficient.
   Coefficient &corrf; ///< Borrowed correction coefficient.

public:
   /**
    * @brief Borrow coefficients for abs(basec_) * corrc_.
    * @param basec_ Coefficient whose absolute value is used.
    * @param corrc_ Multiplicative correction coefficient.
    */
   PProductCoefficient(Coefficient& basec_, Coefficient& corrc_)
      : basef(basec_), corrf(corrc_) { }

   /**
    * @brief Evaluate the scalar coefficient at an integration point.
    * @param T Element transformation.
    * @param ip Reference integration point.
    * @return Coefficient value.
    */
   real_t Eval(ElementTransformation &T,
               const IntegrationPoint &ip) override
   {
      T.SetIntPoint(&ip);
      real_t u = basef.Eval(T,ip);
      real_t c = corrf.Eval(T,ip);
      if (u<0.0) { u*=-1.0; }
      return u*c;
   }
};


// Formulation for the ScreenedPoisson equation. The positive part of the input
// coefficient supply unit volumetric loading, the negative part - negative unit
// volumetric loading. The parameter rh is the radius of a linear cone filter
// which will deliver similar smoothing effect as the Screened Poisson
// equation. It determines the length scale of the smoothing.
/** @brief Screened-Poisson energy, residual, and Jacobian integrator. */
class ScreenedPoisson: public NonlinearFormIntegrator
{
protected:
   real_t diffcoef; ///< Squared screened-Poisson smoothing length.
   Coefficient *func; ///< Input coefficient; ownership depends on the integrator.

public:
   /**
    * @brief Configure screened-Poisson forcing and smoothing scale.
    * @param nfunc Borrowed source coefficient.
    * @param rh Equivalent cone-filter radius.
    */
   ScreenedPoisson(Coefficient &nfunc, real_t rh):func(&nfunc)
   {
      real_t rd=rh/(2*std::sqrt(3.0));
      diffcoef= rd*rd;
   }

   /**
    * @brief Destroy the integrator without deleting the borrowed source.
    */
   ~ScreenedPoisson() { }

   /**
    * @brief Replace the borrowed source coefficient.
    * @param nfunc New input coefficient.
    */
   void SetInput(Coefficient &nfunc) { func = &nfunc; }

   /**
    * @brief Evaluate the element contribution to the energy.
    * @param el Finite element.
    * @param trans Element transformation.
    * @param elfun Element solution values.
    * @return Element energy.
    */
   real_t GetElementEnergy(const FiniteElement &el,
                           ElementTransformation &trans,
                           const Vector &elfun) override;

   /**
    * @brief Assemble the element residual.
    * @param el Finite element.
    * @param trans Element transformation.
    * @param elfun Element solution values.
    * @param[out] elvect Element residual vector.
    */
   void AssembleElementVector(const FiniteElement &el,
                              ElementTransformation &trans,
                              const Vector &elfun,
                              Vector &elvect) override;

   /**
    * @brief Assemble the element Jacobian.
    * @param el Finite element.
    * @param trans Element transformation.
    * @param elfun Element solution values.
    * @param[out] elmat Element Jacobian matrix.
    */
   void AssembleElementGrad(const FiniteElement &el,
                            ElementTransformation &trans,
                            const Vector &elfun,
                            DenseMatrix &elmat) override;
};


/** @brief Regularized p-Laplacian energy, residual, and Jacobian integrator. */
class PUMPLaplacian: public NonlinearFormIntegrator
{

protected:
   Coefficient *func; ///< Input coefficient; ownership depends on the integrator.
   VectorCoefficient *fgrad; ///< Input gradient coefficient.
   bool ownership; ///< Whether to delete input coefficients on destruction.
   real_t pp; ///< p-Laplacian power.
   real_t ee; ///< Regularization parameter.

public:
   /// The VectorCoefficent should contain a vector with entries:
   /// [0] - derivative with respect to x
   /// [1] - derivative with respect to y
   /// [2] - derivative with respect to z
   /**
    * @brief Configure source, gradient, and coefficient ownership.
    * @param nfunc Source coefficient.
    * @param nfgrad Gradient coefficient with one component per spatial dimension.
    * @param ownership_ Whether the destructor deletes both coefficients.
    */
   PUMPLaplacian(Coefficient *nfunc, VectorCoefficient *nfgrad,
                 bool ownership_=true)
      : func(nfunc), fgrad(nfgrad), ownership(ownership_), pp(2.0), ee(1e-7) { }

   /**
    * @brief Set the p-Laplacian exponent.
    * @param pp_ Exponent.
    */
   void SetPower(real_t pp_) { pp = pp_; }
   /**
    * @brief Set the regularization parameter.
    * @param ee_ Regularization value.
    */
   void SetReg(real_t ee_)   { ee = ee_; }

   /**
    * @brief Delete the source and gradient coefficients when ownership is enabled.
    */
   virtual ~PUMPLaplacian()
   {
      if (ownership)
      {
         delete func;
         delete fgrad;
      }
   }

   /**
    * @brief Evaluate the element contribution to the energy.
    * @param el Finite element.
    * @param trans Element transformation.
    * @param elfun Element solution values.
    * @return Element energy.
    */
   real_t GetElementEnergy(const FiniteElement &el,
                           ElementTransformation &trans,
                           const Vector &elfun) override;

   /**
    * @brief Assemble the element residual.
    * @param el Finite element.
    * @param trans Element transformation.
    * @param elfun Element solution values.
    * @param[out] elvect Element residual vector.
    */
   void AssembleElementVector(const FiniteElement &el,
                              ElementTransformation &trans,
                              const Vector &elfun,
                              Vector &elvect) override;

   /**
    * @brief Assemble the element Jacobian.
    * @param el Finite element.
    * @param trans Element transformation.
    * @param elfun Element solution values.
    * @param[out] elmat Element Jacobian matrix.
    */
   void AssembleElementGrad(const FiniteElement &el,
                            ElementTransformation &trans,
                            const Vector &elfun,
                            DenseMatrix &elmat) override;
};

// Low-pass filter based on the Screened Poisson equation.
// B. S. Lazarov, O. Sigmund: "Filters in topology optimization based on
// Helmholtz-type differential equations", DOI:10.1002/nme.3072.
/** @brief Low-pass filter based on a screened-Poisson solve. */
class PDEFilter
{
public:
   /**
    * @brief Construct a parallel screened-Poisson filter and linear solver.
    * @param mesh Mesh defining the filter domain; must outlive the filter.
    * @param rh Equivalent cone-filter radius.
    * @param order H1 polynomial order.
    * @param maxiter Maximum GMRES iterations.
    * @param rtol Relative solver tolerance.
    * @param atol Absolute solver tolerance.
    * @param print_lv Solver and preconditioner output level.
    */
   PDEFilter(ParMesh &mesh, real_t rh, int order = 2,
             int maxiter = 100, real_t rtol = 1e-12,
             real_t atol = 1e-15, int print_lv = 0)
      : rr(rh),
        fecp(order, mesh.Dimension()),
        fesp(&mesh, &fecp, 1),
        gf(&fesp)
   {
      sv = fesp.NewTrueDofVector();

      nf = new ParNonlinearForm(&fesp);
      prec = new HypreBoomerAMG();
      prec->SetPrintLevel(print_lv);

      gmres = new GMRESSolver(mesh.GetComm());

      gmres->SetAbsTol(atol);
      gmres->SetRelTol(rtol);
      gmres->SetMaxIter(maxiter);
      gmres->SetPrintLevel(print_lv);
      gmres->SetPreconditioner(*prec);

      sint=nullptr;
   }

   /**
    * @brief Release the owned solver, preconditioner, nonlinear form, and vector.
    */
   ~PDEFilter()
   {
      delete gmres;
      delete prec;
      delete nf;
      delete sv;
   }

   /**
    * @brief Filter a grid function through a coefficient wrapper.
    * @param func Input field.
    * @param[out] ffield Filtered output on its initialized space.
    */
   void Filter(ParGridFunction &func, ParGridFunction &ffield)
   {
      GridFunctionCoefficient gfc(&func);
      Filter(gfc, ffield);
   }

   /**
    * @brief Solve the screened-Poisson filtering problem.
    * @param func Input coefficient.
    * @param[out] ffield Filtered output on its initialized space.
    */
   void Filter(Coefficient &func, ParGridFunction &ffield);

private:
   const real_t rr; ///< Equivalent cone-filter radius.
   H1_FECollection fecp; ///< Owned scalar H1 finite-element collection.
   ParFiniteElementSpace fesp; ///< Filter space on the borrowed mesh.
   ParGridFunction gf; ///< Filter solution workspace.

   ParNonlinearForm* nf; ///< Owned screened-Poisson nonlinear form.
   HypreBoomerAMG* prec; ///< Owned AMG preconditioner.
   GMRESSolver *gmres; ///< Owned linear solver.
   HypreParVector *sv; ///< Owned true-DOF solution vector.

   ScreenedPoisson* sint; ///< Integrator owned by the nonlinear form.
};

} // namespace common

} // namespace mfem

#endif // MFEM_USE_MPI
#endif
