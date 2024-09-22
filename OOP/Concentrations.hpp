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
    void InitializeCnP(ParFiniteElementSpace *fespace);
    void InitializeCnE(ParFiniteElementSpace *fespace);

    void TimeStepCnP(ParFiniteElementSpace *fespace);
    void TimeStepCnE(ParFiniteElementSpace *fespace);

    void SetupBoundaryConditions(ParFiniteElementSpace *fespace);

    // void TestFESpace(ParFiniteElementSpace *fespace);

    mfem::ParGridFunction PeR;
    mfem::GridFunctionCoefficient matCoef_R;

    // std::shared_ptr<ParFiniteElementSpace> fespace;
    // void TestFESpace(std::shared_ptr<ParFiniteElementSpace> fespace);


private:

    MeshHandler &mesh_handler;
    Reaction reaction;

    // Internal methods for different operations
    void CreateCnE(mfem::ParGridFunction &Cn, double initial_value);
    void Lithiation(mfem::ParGridFunction &Cn, double initial_value, ParFiniteElementSpace *fespace);
    void SBM_Matrix(mfem::ParGridFunction &psx, HypreParMatrix &Mmat, ParFiniteElementSpace *fespace);
    void Solver(HypreParMatrix &Mmat, CGSolver &M_solver);
    void ImposeNeumannBC(mfem::ParGridFunction &PGF, mfem::ParGridFunction &psx);
    void SetupRx(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, double value, GridFunctionCoefficient cAx);
    void ForceTerm(ParFiniteElementSpace *fespace, GridFunctionCoefficient cXx, mfem::ParLinearForm &Fxx, Array<int> boundary, ConstantCoefficient m);
    void TotalReaction(mfem::ParGridFunction &Rx, double xCrnt);
    void EnsureValidBoundaryAndFESpace();
    void DebugBoundaryArray(const Array<int> &boundary);

    // Member variables
    mfem::ParFiniteElementSpace* fespace;           // Finite element space
    mfem::ParGridFunction psi;                     // Psi grid function from MeshHandler
    mfem::ParGridFunction pse;                     // Pse grid function from MeshHandler
    mfem::ParGridFunction TmpF;
    mfem::ParMesh* pmesh;

    const mfem::Vector& EVol;                       // Element volumes from MeshHandler
    double gtPsi;                                   // Total Psi from MeshHandler
    double infx;

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
    mfem::ParLinearForm Fet;
    // mfem::ParGridFunction PeR;

    // GridFunctionCoefficient matCoef_R;
    GridFunctionCoefficient cAe;

    std::unique_ptr<mfem::ProductCoefficient> m_nbcCoef;    
    mfem::Array<int> nbc_w_bdr;       

    double eCrnt;

};

#endif // CONCENTRATIONS_HPP
