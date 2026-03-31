// // #ifndef CNE_HPP
// // #define CNE_HPP

// // #include "Concentrations_Base.hpp"
// // #include "SimTypes.hpp"

// // class Initialize_Geometry;
// // class Domain_Parameters;
// // class BoundaryConditions;

// // /**
// //  * @class CnE
// //  * @brief Electrolyte concentration solver using diffusion + reaction coupling.
// //  *
// //  * Implements a Crank–Nicolson (CN) update for the electrolyte concentration
// //  * field. The class assembles mass and stiffness operators, applies Neumann
// //  * boundary fluxes, handles reaction source terms, and performs time-stepping
// //  * using MFEM/Hypre parallel solvers.
// //  *
// //  * Supports both full-cell (anode + cathode) and half-cell configurations.
// //  */
// // class CnE : public ConcentrationBase {
// // public:

// //     /**
// //      * @brief Construct the electrolyte concentration solver.
// //      *
// //      * @param geo  Geometry handler (mesh, FE spaces, DOF mappings).
// //      * @param para Domain parameters (material, time-step constants).
// //      * @param bc   Boundary conditions handler (Neumann/Dirichlet markers).
// //      * @param mode Cell mode (full cell, anode-only, cathode-only).
// //      */
// //     CnE(Initialize_Geometry &geo,
// //         Domain_Parameters  &para,
// //         BoundaryConditions &bc,
// //         sim::CellMode mode);

// //     /**
// //      * @brief Initialize electrolyte concentration and assemble operators.
// //      *
// //      * Performs:
// //      * - Initial concentration assignment  
// //      * - Assembly of mass and stiffness matrices  
// //      * - Setup of boundary forcing terms (Neumann BCs)  
// //      * - Initialization of CG solver and preconditioner  
// //      *
// //      * @param Cn            Electrolyte concentration field (in/out).
// //      * @param initial_value Initial scalar concentration value.
// //      * @param psx           Phase-field ψ for electrolyte-region masking.
// //      */
// //     void SetupField(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx);

// //     /**
// //      * @brief Advance electrolyte concentration by one timestep for half cell (single Rx source).
// //      *
// //      * Constructs CN left/right operators, assembles updated reaction forcing,
// //      * applies Neumann BCs, solves for the new concentration vector, and restores
// //      * the field to \p Cn (ParGridFunction).
// //      *
// //      * @param Rx  Reaction field (combined or single-source).
// //      * @param Cn  Electrolyte concentration field (input/output).
// //      * @param psx Phase-field ψ for diffusivity/BC masking.
// //      */
// //     void UpdateConcentration(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);

// //     /**
// //      * @brief Advance electrolyte concentration with separate anode/cathode sources for full cell.
// //      *
// //      * Same CN update as the single-source version, but allows separate reaction
// //      * contributions from the anode and cathode before combining them.
// //      *
// //      * @param RxC Reaction field from the cathode.
// //      * @param RxA Reaction field from the anode.
// //      * @param Cn  Electrolyte concentration field (input/output).
// //      * @param psx Phase-field ψ for diffusivity/BC masking.
// //      */
// //     void UpdateConcentration(mfem::ParGridFunction &RxC, mfem::ParGridFunction &RxA,
// //                              mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);

// //     /**
// //      * @brief Enforce salt conservation.
// //      *
// //      * Applies a correction ensuring the electrolyte salt inventory remains
// //      * conserved over the domain.
// //      *
// //      * @param Cn  Electrolyte concentration field to normalize.
// //      * @param psx Phase-field ψ mask.
// //      */
// //     void SaltConservation(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);

// //     // -------------------------------------------------------------------------
// //     // Geometry / spaces
// //     // -------------------------------------------------------------------------
// //     Initialize_Geometry &geometry;            ///< Geometry & mesh container.
// //     Domain_Parameters   &domain_parameters;   ///< Electrolyte material parameters.
// //     BoundaryConditions  &boundary_conditions; ///< Boundary-condition handler.
// //     sim::CellMode        mode_;               ///< Cell mode (full/anode/cathode).

// //     std::shared_ptr<mfem::ParFiniteElementSpace> fespace; ///< Parallel FE space.

// //     FEMOperators fem; ///< Operator construction helpers (mass, stiffness).
// //     Utils        utils; ///< Utility routines (clamping, checks, etc.).

// //     mfem::HypreParVector CeVn; ///< Concentration at next time-step (true DoFs).

// //     // -------------------------------------------------------------------------
// //     // Boundary forcing
// //     // -------------------------------------------------------------------------
// //     mfem::ParLinearForm Fet;                 ///< Domain and boundary forcing form.
// //     mfem::Array<int>    nbc_bdr;             ///< Neumann boundary markers.
// //     std::unique_ptr<mfem::ProductCoefficient> m_nbcCoef; ///< Multiply coefficient for BCs.
// //     mfem::ConstantCoefficient nbcCoef;       ///< Scalar Neumann coefficient.
// //     mfem::GridFunctionCoefficient matCoef_R; ///< Coefficient wrapper for PeR.
// //     mfem::Array<int> boundary_dofs;          ///< Boundary true DOFs.

// // private:

// //     // -------------------------------------------------------------------------
// //     // Operators (mass, stiffness, forcing)
// //     // -------------------------------------------------------------------------
// //     std::unique_ptr<mfem::ParBilinearForm> Me_init; ///< Mass bilinear form.
// //     std::unique_ptr<mfem::ParLinearForm>   Be_init; ///< Domain/boundary forcing form.

// //     mfem::HypreParMatrix Mmate; ///< Mass matrix M.
// //     mfem::HypreParMatrix Kmate; ///< Stiffness matrix K (diffusion).
// //     std::unique_ptr<mfem::ParBilinearForm> Ke2; ///< Diffusion stiffness form.

// //     // -------------------------------------------------------------------------
// //     // Solver / preconditioner
// //     // -------------------------------------------------------------------------
// //     mfem::CGSolver      Me_solver; ///< CG solver for M x = b.
// //     mfem::HypreSmoother Me_prec;   ///< Hypre smoother / preconditioner.

// //     // -------------------------------------------------------------------------
// //     // Coefficients / fields
// //     // -------------------------------------------------------------------------
// //     mfem::GridFunctionCoefficient cAe; ///< Reaction coefficient wrapper.
// //     mfem::GridFunctionCoefficient cDe; ///< Diffusivity coefficient wrapper.

// //     mfem::ParGridFunction De;  ///< Diffusivity field for electrolyte.
// //     mfem::ParGridFunction Rxe; ///< Reaction source field (electrolyte).
// //     mfem::ParGridFunction PeR; ///< Reaction potential field.

// //     // -------------------------------------------------------------------------
// //     // State values (diagnostics)
// //     // -------------------------------------------------------------------------
// //     double eCrnt = 0.0; ///< Reaction current.
// //     double infx  = 0.0; ///< Boundary flux value.
// //     double L_w   = 0.0; ///< Characteristic electrolyte width.

// //     // -------------------------------------------------------------------------
// //     // Vectors (true DoF storage)
// //     // -------------------------------------------------------------------------
// //     mfem::HypreParVector Feb; ///< Assembled RHS vector.
// //     mfem::HypreParVector X1v; ///< Workspace vector.

// //     std::unique_ptr<mfem::HypreParMatrix> TmatR; ///< CN right operator  (M - 0.5 dt K).
// //     std::unique_ptr<mfem::HypreParMatrix> TmatL; ///< CN left operator   (M + 0.5 dt K).

// //     mfem::HypreParVector CeV0; ///< Concentration at current step.
// //     mfem::HypreParVector RHCe; ///< RHS after CN assembly (TmatR + forcing).
// //     mfem::HypreParVector PsVc; ///< ψ-field mask in true-DoF representation.

// // };

// // #endif // CNE_HPP

// #ifndef CNE_HPP
// #define CNE_HPP

// #include "Concentrations_Base.hpp"
// #include "SimTypes.hpp"

// class Initialize_Geometry;
// class Domain_Parameters;
// class BoundaryConditions;

// /**
//  * @class CnE
//  * @brief Electrolyte concentration solver using diffusion + reaction coupling.
//  *
//  * Implements a Crank–Nicolson (CN) update for the electrolyte concentration
//  * field. The class assembles mass and stiffness operators, applies Neumann
//  * boundary fluxes, handles reaction source terms, and performs time-stepping
//  * using MFEM/Hypre parallel solvers.
//  *
//  * Supports both full-cell (anode + cathode) and half-cell configurations.
//  */
// class CnE : public ConcentrationBase {
// public:

//     /**
//      * @brief Construct the electrolyte concentration solver.
//      *
//      * @param geo  Geometry handler (mesh, FE spaces, DOF mappings).
//      * @param para Domain parameters (material, time-step constants).
//      * @param bc   Boundary conditions handler (Neumann/Dirichlet markers).
//      * @param mode Cell mode (full cell, anode-only, cathode-only).
//      */
//     CnE(Initialize_Geometry &geo,
//         Domain_Parameters  &para,
//         BoundaryConditions &bc,
//         sim::CellMode mode);

//     /**
//      * @brief Initialize electrolyte concentration and assemble operators.
//      *
//      * Performs:
//      * - Initial concentration assignment  
//      * - Assembly of mass and stiffness matrices  
//      * - Setup of boundary forcing terms (Neumann BCs)  
//      * - Initialization of CG solver and preconditioner  
//      *
//      * @param Cn            Electrolyte concentration field (in/out).
//      * @param initial_value Initial scalar concentration value.
//      * @param psx           Phase-field ψ for electrolyte-region masking.
//      * @param gtPsx         Global normalization of ψ (used for boundary conditions).
//      */
//     void SetupField(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx, double gtPsx);
//     /**
//      * @brief Advance electrolyte concentration by one timestep for half cell (single Rx source).
//      *
//      * Constructs CN left/right operators, assembles updated reaction forcing,
//      * applies Neumann BCs, solves for the new concentration vector, and restores
//      * the field to \p Cn (ParGridFunction).
//      *
//      * @param Rx  Reaction field (combined or single-source).
//      * @param Cn  Electrolyte concentration field (input/output).
//      * @param psx Phase-field ψ for diffusivity/BC masking.
//      * @param gtPsx Global normalization of ψ (used for boundary conditions).
//      */
//     void UpdateConcentration(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, double gtPsx);

//     /**
//      * @brief Advance electrolyte concentration with separate anode/cathode sources for full cell.
//      *
//      * Same CN update as the single-source version, but allows separate reaction
//      * contributions from the anode and cathode before combining them.
//      *
//      * @param RxC Reaction field from the cathode.
//      * @param RxA Reaction field from the anode.
//      * @param Cn  Electrolyte concentration field (input/output).
//      * @param psx Phase-field ψ for diffusivity/BC masking.
//      */
//     void UpdateConcentration(mfem::ParGridFunction &RxC, mfem::ParGridFunction &RxA,
//                              mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);


//     /**
//      * @brief Advance electrolyte concentration with particle interactions.
//      *
//      * Updates the electrolyte concentration field using a Crank–Nicolson scheme
//      * that includes particle interaction effects (e.g., AB, AC).
//      *
//      * @param Rx  Reaction field (combined or single-source).
//      * @param Cn  Electrolyte concentration field (input/output).
//      * @param psx Phase-field ψ for diffusivity/BC masking.
//      * @param gtPsx Global normalization of ψ (used for boundary conditions).
//      * @param weight_elec Weighting function for electrode contributions.
//      * @param sum_AB Summation of particle interactions for AB.
//      * @param weight_AB Weighting function for AB interactions.
//      * @param grad_AB Gradient of AB interactions.
//      * @param sum_AC Summation of particle interactions for AC.
//      * @param weight_AC Weighting function for AC interactions.
//      * @param grad_AC Gradient of AC interactions.
//      * @param mu_A Chemical potential of particle A.
//      * @param mu_B Chemical potential of particle B.
//      * @param mu_C Chemical potential of particle C.
//      * @param mu_D Chemical potential of particle D.
//      */
//     void UpdateConcentration(mfem::ParGridFunction &Rx,
//                              mfem::ParGridFunction &Cn,
//                              mfem::ParGridFunction &psx,
//                              double gtPsx,
//                              mfem::ParGridFunction &weight_elec,
//                              mfem::ParGridFunction &sum_AB,
//                              mfem::ParGridFunction &weight_AB,
//                              mfem::ParGridFunction &grad_AB,
//                              mfem::ParGridFunction &sum_AC,
//                              mfem::ParGridFunction &weight_AC,
//                              mfem::ParGridFunction &grad_AC,
//                              mfem::ParGridFunction &mu_A, mfem::ParGridFunction &mu_B, mfem::ParGridFunction &mu_C, mfem::ParGridFunction &mu_D, mfem::ParGridFunction &psiA, mfem::ParGridFunction &psiB, mfem::ParGridFunction &psiC);

//     /**
//      * @brief Enforce salt conservation.
//      *
//      * Applies a correction ensuring the electrolyte salt inventory remains
//      * conserved over the domain.
//      *
//      * @param Cn  Electrolyte concentration field to normalize.
//      * @param psx Phase-field ψ mask.
//      */
//     void SaltConservation(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);

//     // -------------------------------------------------------------------------
//     // Geometry / spaces
//     // -------------------------------------------------------------------------
//     Initialize_Geometry &geometry;            ///< Geometry & mesh container.
//     Domain_Parameters   &domain_parameters;   ///< Electrolyte material parameters.
//     BoundaryConditions  &boundary_conditions; ///< Boundary-condition handler.
//     sim::CellMode        mode_;               ///< Cell mode (full/anode/cathode).

//     std::shared_ptr<mfem::ParFiniteElementSpace> fespace; ///< Parallel FE space.

//     FEMOperators fem; ///< Operator construction helpers (mass, stiffness).
//     Utils        utils; ///< Utility routines (clamping, checks, etc.).

//     mfem::HypreParVector CeVn; ///< Concentration at next time-step (true DoFs).

//     // -------------------------------------------------------------------------
//     // Boundary forcing
//     // -------------------------------------------------------------------------
//     mfem::ParLinearForm Fet;                 ///< Domain and boundary forcing form.
//     mfem::Array<int>    nbc_bdr;             ///< Neumann boundary markers.
//     std::unique_ptr<mfem::ProductCoefficient> m_nbcCoef; ///< Multiply coefficient for BCs.
//     mfem::ConstantCoefficient nbcCoef;       ///< Scalar Neumann coefficient.
//     mfem::GridFunctionCoefficient matCoef_R; ///< Coefficient wrapper for PeR.
//     mfem::Array<int> boundary_dofs;          ///< Boundary true DOFs.

// private:

//     // -------------------------------------------------------------------------
//     // Operators (mass, stiffness, forcing)
//     // -------------------------------------------------------------------------
//     std::unique_ptr<mfem::ParBilinearForm> Me_init; ///< Mass bilinear form.
//     std::unique_ptr<mfem::ParLinearForm>   Be_init; ///< Domain/boundary forcing form.

//     mfem::HypreParMatrix Mmate; ///< Mass matrix M.
//     mfem::HypreParMatrix Kmate; ///< Stiffness matrix K (diffusion).
//     std::unique_ptr<mfem::ParBilinearForm> Ke2; ///< Diffusion stiffness form.

//     // -------------------------------------------------------------------------
//     // Solver / preconditioner
//     // -------------------------------------------------------------------------
//     mfem::CGSolver      Me_solver; ///< CG solver for M x = b.
//     mfem::HypreSmoother Me_prec;   ///< Hypre smoother / preconditioner.

//     // -------------------------------------------------------------------------
//     // Coefficients / fields
//     // -------------------------------------------------------------------------
//     mfem::GridFunctionCoefficient cAe; ///< Reaction coefficient wrapper.
//     mfem::GridFunctionCoefficient cDe; ///< Diffusivity coefficient wrapper.

//     mfem::ParGridFunction De;  ///< Diffusivity field for electrolyte.
//     mfem::ParGridFunction Rxe; ///< Reaction source field (electrolyte).
//     mfem::ParGridFunction PeR; ///< Reaction potential field.

//     // -------------------------------------------------------------------------
//     // State values (diagnostics)
//     // -------------------------------------------------------------------------
//     double eCrnt = 0.0; ///< Reaction current.
//     double infx  = 0.0; ///< Boundary flux value.
//     double L_w   = 0.0; ///< Characteristic electrolyte width.

//     // -------------------------------------------------------------------------
//     // Vectors (true DoF storage)
//     // -------------------------------------------------------------------------
//     mfem::HypreParVector Feb; ///< Assembled RHS vector.
//     mfem::HypreParVector X1v; ///< Workspace vector.

//     std::unique_ptr<mfem::HypreParMatrix> TmatR; ///< CN right operator  (M - 0.5 dt K).
//     std::unique_ptr<mfem::HypreParMatrix> TmatL; ///< CN left operator   (M + 0.5 dt K).

//     mfem::HypreParVector CeV0; ///< Concentration at current step.
//     mfem::HypreParVector RHCe; ///< RHS after CN assembly (TmatR + forcing).
//     mfem::HypreParVector PsVc; ///< ψ-field mask in true-DoF representation.

// };

// #endif // CNE_HPP

#ifndef CNE_HPP
#define CNE_HPP

#include "Concentrations_Base.hpp"
#include "SimTypes.hpp"

class Initialize_Geometry;
class Domain_Parameters;
class BoundaryConditions;

/**
 * @class CnE
 * @brief Electrolyte concentration solver using diffusion + reaction coupling.
 *
 * Implements a Crank–Nicolson (CN) update for the electrolyte concentration
 * field. The class assembles mass and stiffness operators, applies Neumann
 * boundary fluxes, handles reaction source terms, and performs time-stepping
 * using MFEM/Hypre parallel solvers.
 *
 * Supports both full-cell (anode + cathode) and half-cell configurations.
 */
class CnE : public ConcentrationBase {
public:

    /**
     * @brief Construct the electrolyte concentration solver.
     *
     * @param geo  Geometry handler (mesh, FE spaces, DOF mappings).
     * @param para Domain parameters (material, time-step constants).
     * @param bc   Boundary conditions handler (Neumann/Dirichlet markers).
     * @param mode Cell mode (full cell, anode-only, cathode-only).
     */
    CnE(Initialize_Geometry &geo,
        Domain_Parameters  &para,
        BoundaryConditions &bc,
        sim::CellMode mode);

    /**
     * @brief Initialize electrolyte concentration and assemble operators.
     *
     * Performs:
     * - Initial concentration assignment  
     * - Assembly of mass and stiffness matrices  
     * - Setup of boundary forcing terms (Neumann BCs)  
     * - Initialization of CG solver and preconditioner  
     *
     * @param Cn            Electrolyte concentration field (in/out).
     * @param initial_value Initial scalar concentration value.
     * @param psx           Phase-field ψ for electrolyte-region masking.
     * @param gtPsx         Global normalization of ψ (used for boundary conditions).
     */
    void SetupField(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx, double gtPsx);
    /**
     * @brief Advance electrolyte concentration by one timestep for half cell (single Rx source).
     *
     * Constructs CN left/right operators, assembles updated reaction forcing,
     * applies Neumann BCs, solves for the new concentration vector, and restores
     * the field to \p Cn (ParGridFunction).
     *
     * @param Rx  Reaction field (combined or single-source).
     * @param Cn  Electrolyte concentration field (input/output).
     * @param psx Phase-field ψ for diffusivity/BC masking.
     * @param gtPsx Global normalization of ψ (used for boundary conditions).
     */
    void UpdateConcentration(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, double gtPsx);

    /**
     * @brief Advance electrolyte concentration with separate anode/cathode sources for full cell.
     *
     * Same CN update as the single-source version, but allows separate reaction
     * contributions from the anode and cathode before combining them.
     *
     * @param RxC Reaction field from the cathode.
     * @param RxA Reaction field from the anode.
     * @param Cn  Electrolyte concentration field (input/output).
     * @param psx Phase-field ψ for diffusivity/BC masking.
     */
    void UpdateConcentration(mfem::ParGridFunction &RxC, mfem::ParGridFunction &RxA,
                             mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);


    /**
     * @brief Advance electrolyte concentration with particle interactions.
     *
     * Updates the electrolyte concentration field using a Crank–Nicolson scheme
     * that includes particle interaction effects (e.g., AB, AC).
     *
     * @param Rx  Reaction field (combined or single-source).
     * @param Cn  Electrolyte concentration field (input/output).
     * @param psx Phase-field ψ for diffusivity/BC masking.
     * @param gtPsx Global normalization of ψ (used for boundary conditions).
     * @param weight_elec Weighting function for electrode contributions.
     * @param sum_AB Summation of particle interactions for AB.
     * @param weight_AB Weighting function for AB interactions.
     * @param grad_AB Gradient of AB interactions.
     * @param sum_AC Summation of particle interactions for AC.
     * @param weight_AC Weighting function for AC interactions.
     * @param grad_AC Gradient of AC interactions.
     * @param mu_A Chemical potential of particle A.
     * @param mu_B Chemical potential of particle B.
     * @param mu_C Chemical potential of particle C.
     * @param mu_D Chemical potential of particle D.
     */
    void UpdateConcentration(mfem::ParGridFunction &Rx,
                             mfem::ParGridFunction &Cn,
                             mfem::ParGridFunction &psx,
                             double gtPsx,
                             mfem::ParGridFunction &weight_elec,
                             mfem::ParGridFunction &sum_AB,
                             mfem::ParGridFunction &weight_AB,
                             mfem::ParGridFunction &grad_AB,
                             mfem::ParGridFunction &sum_AC,
                             mfem::ParGridFunction &weight_AC,
                             mfem::ParGridFunction &grad_AC,
                             mfem::ParGridFunction &mu_A, mfem::ParGridFunction &mu_B, mfem::ParGridFunction &mu_C, mfem::ParGridFunction &mu_D, mfem::ParGridFunction &psiA, mfem::ParGridFunction &psiB, mfem::ParGridFunction &psiC);

    /**
     * @brief Enforce salt conservation.
     *
     * Applies a correction ensuring the electrolyte salt inventory remains
     * conserved over the domain.
     *
     * @param Cn  Electrolyte concentration field to normalize.
     * @param psx Phase-field ψ mask.
     */
    void SaltConservation(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);

    // -------------------------------------------------------------------------
    // Geometry / spaces
    // -------------------------------------------------------------------------
    Initialize_Geometry &geometry;            ///< Geometry & mesh container.
    Domain_Parameters   &domain_parameters;   ///< Electrolyte material parameters.
    BoundaryConditions  &boundary_conditions; ///< Boundary-condition handler.
    sim::CellMode        mode_;               ///< Cell mode (full/anode/cathode).

    std::shared_ptr<mfem::ParFiniteElementSpace> fespace; ///< Parallel FE space.

    FEMOperators fem; ///< Operator construction helpers (mass, stiffness).
    Utils        utils; ///< Utility routines (clamping, checks, etc.).

    mfem::HypreParVector CeVn; ///< Concentration at next time-step (true DoFs).

    // -------------------------------------------------------------------------
    // Boundary forcing
    // -------------------------------------------------------------------------
    mfem::ParLinearForm Fet;                 ///< Domain and boundary forcing form.
    mfem::Array<int>    nbc_bdr;             ///< Neumann boundary markers.
    std::unique_ptr<mfem::ProductCoefficient> m_nbcCoef; ///< Multiply coefficient for BCs.
    mfem::ConstantCoefficient nbcCoef;       ///< Scalar Neumann coefficient.
    mfem::GridFunctionCoefficient matCoef_R; ///< Coefficient wrapper for PeR.
    mfem::Array<int> boundary_dofs;          ///< Boundary true DOFs.

private:

    // -------------------------------------------------------------------------
    // Operators (mass, stiffness, forcing)
    // -------------------------------------------------------------------------
    std::unique_ptr<mfem::ParBilinearForm> Me_init; ///< Mass bilinear form.
    std::unique_ptr<mfem::ParLinearForm>   Be_init; ///< Domain/boundary forcing form.

    mfem::HypreParMatrix Mmate; ///< Mass matrix M.
    mfem::HypreParMatrix Kmate; ///< Stiffness matrix K (diffusion).
    std::unique_ptr<mfem::ParBilinearForm> Ke2; ///< Diffusion stiffness form.

    // -------------------------------------------------------------------------
    // Solver / preconditioner
    // -------------------------------------------------------------------------
    mfem::CGSolver      Me_solver; ///< CG solver for M x = b.
    mfem::HypreSmoother Me_prec;   ///< Hypre smoother / preconditioner.

    // -------------------------------------------------------------------------
    // Coefficients / fields
    // -------------------------------------------------------------------------
    mfem::GridFunctionCoefficient cAe; ///< Reaction coefficient wrapper.
    mfem::GridFunctionCoefficient cDe; ///< Diffusivity coefficient wrapper.

    mfem::ParGridFunction De;  ///< Diffusivity field for electrolyte.
    mfem::ParGridFunction Rxe; ///< Reaction source field (electrolyte).
    mfem::ParGridFunction PeR; ///< Reaction potential field.

    // -------------------------------------------------------------------------
    // State values (diagnostics)
    // -------------------------------------------------------------------------
    double eCrnt = 0.0; ///< Reaction current.
    double infx  = 0.0; ///< Boundary flux value.
    double L_w   = 0.0; ///< Characteristic electrolyte width.

    // -------------------------------------------------------------------------
    // Vectors (true DoF storage)
    // -------------------------------------------------------------------------
    mfem::HypreParVector Feb; ///< Assembled RHS vector.
    mfem::HypreParVector X1v; ///< Workspace vector.

    std::unique_ptr<mfem::HypreParMatrix> TmatR; ///< CN right operator  (M - 0.5 dt K).
    std::unique_ptr<mfem::HypreParMatrix> TmatL; ///< CN left operator   (M + 0.5 dt K).

    mfem::HypreParVector CeV0; ///< Concentration at current step.
    mfem::HypreParVector RHCe; ///< RHS after CN assembly (TmatR + forcing).
    mfem::HypreParVector PsVc; ///< ψ-field mask in true-DoF representation.

};

#endif // CNE_HPP