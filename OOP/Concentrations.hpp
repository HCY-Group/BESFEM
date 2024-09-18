#ifndef CONCENTRATIONS_HPP
#define CONCENTRATIONS_HPP

#include "mfem.hpp"
#include "Mesh_Handler.hpp"
#include "Constants.hpp"
#include "Reaction.hpp"
#include <memory>

class Concentrations {
public:
    // Constructor that takes a MeshHandler reference
    Concentrations(MeshHandler &mesh_handler);

    // Initialization method
    void InitializeCnP();
    void InitializeCnE();

    void TimeStepCnP();
    void TimeStepCnE();

private:

    MeshHandler &mesh_handler;
    Reaction reaction;

    // Internal methods for different operations
    void CreateCnE(mfem::ParGridFunction &Cn, double initial_value);
    void Lithiation(mfem::ParGridFunction &Cn, double initial_value);
    void SBM_Matrix(mfem::ParGridFunction &psx, HypreParMatrix &Mmat);
    void Solver(HypreParMatrix &Mmat, CGSolver &M_solver);
    void SetupBoundaryConditions();
    void SetupRx(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, double value, GridFunctionCoefficient cAx);
    void ForceTerm(GridFunctionCoefficient cXx, mfem::ParLinearForm &Fxx);
    void TotalReaction(mfem::ParGridFunction &Rx, double xCrnt);

    // Member variables
    mfem::ParFiniteElementSpace* fespace;           // Finite element space
    mfem::ParGridFunction psi;                     // Psi grid function from MeshHandler
    mfem::ParGridFunction pse;                     // Pse grid function from MeshHandler
    mfem::ParGridFunction TmpF;

    const mfem::Vector& EVol;                       // Element volumes from MeshHandler
    double gtPsi;                                   // Total Psi from MeshHandler

    int nE;                                         // Number of elements
    int nC;                                         // Number of corners (assuming this is number of corners)
    int nV;                                         // Number of vertices

    Array<int> boundary_dofs;

    CGSolver Mp_solver;
    CGSolver Me_solver;

    std::unique_ptr<mfem::ParGridFunction> CnP;     // Concentration CnP (ParGridFunction)
    std::unique_ptr<mfem::ParGridFunction> CnE;     // Concentration CnE (ParGridFunction)

    HypreParMatrix Mmatp;
    HypreParMatrix Mmate;

    std::unique_ptr<mfem::ParGridFunction> Rxc;     
    std::unique_ptr<mfem::ParGridFunction> Rxe;     

    std::unique_ptr<mfem::ParLinearForm> Bc2;
    std::unique_ptr<mfem::ParLinearForm> Be2;

    mfem::ParLinearForm Fct;

    double eCrnt;

};

#endif // CONCENTRATIONS_HPP
