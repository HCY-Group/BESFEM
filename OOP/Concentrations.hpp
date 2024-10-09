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
    void InitializeCnP(std::shared_ptr<ParFiniteElementSpace> fespace);
    void InitializeCnE(std::shared_ptr<ParFiniteElementSpace> fespace);

    void TimeStepCnP(std::shared_ptr<ParFiniteElementSpace> fespace);
    void TimeStepCnE(std::shared_ptr<ParFiniteElementSpace> fespace);

    void SetupBoundaryConditions(std::shared_ptr<ParFiniteElementSpace> fespace);

    // void TestFESpace(ParFiniteElementSpace *fespace);

    mfem::ParGridFunction PeR;
    mfem::GridFunctionCoefficient matCoef_R;

    Array<int> nbc_w_bdr;

    

    mfem::ParFiniteElementSpace* GetFESpace() const { return fespace.get(); }

    // std::shared_ptr<ParFiniteElementSpace> fespace;
    // void TestFESpace(std::shared_ptr<ParFiniteElementSpace> fespace);


private:

    MeshHandler &mesh_handler;
    Reaction reaction;

    // Internal methods for different operations
    void CreateCnE(mfem::ParGridFunction &Cn, double initial_value);
    void Lithiation(mfem::ParGridFunction &Cn, double initial_value, std::shared_ptr<ParFiniteElementSpace> fespace);
    void LithiationCalculation(mfem::ParGridFunction &Cn, std::shared_ptr<ParFiniteElementSpace> fespace);
    void SBM_Matrix(mfem::ParGridFunction &psx, std::shared_ptr<HypreParMatrix> &Mmat, std::shared_ptr<ParFiniteElementSpace> fespace);
    // void Solver(HypreParMatrix &Mmat, CGSolver &M_solver);
    void Solver(std::shared_ptr<HypreParMatrix> &Mmat, std::shared_ptr<CGSolver> &solver);
    void ImposeNeumannBC(mfem::ParGridFunction &PGF, mfem::ParGridFunction &psx);
    void SetupRx(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, double value, GridFunctionCoefficient cAx);
    void ForceTerm(std::shared_ptr<ParFiniteElementSpace> fespace, GridFunctionCoefficient cXx, mfem::ParLinearForm &Fxx, Array<int> boundary, ConstantCoefficient m, bool apply_boundary_conditions);
    void TotalReaction(mfem::ParGridFunction &Rx, double xCrnt);
    void EnsureValidBoundaryAndFESpace();
    void DebugBoundaryArray(const Array<int> &boundary);
    std::shared_ptr<mfem::GridFunctionCoefficient> Diffusivity(mfem::ParGridFunction &psx, mfem::ParGridFunction &Cn, bool particle_electrolyte );
    void K_Matrix(Array<int> boundary, mfem::ParGridFunction &Cn, ParLinearForm &Fxx, std::shared_ptr<HypreParMatrix> &Kmatx, HypreParVector &X1v, HypreParVector &Fxb, GridFunctionCoefficient *cDx);
    
    // Member variables
    // mfem::ParFiniteElementSpace* fespace;           // Finite element space

    std::shared_ptr<mfem::ParFiniteElementSpace> fespace;

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

    std::shared_ptr<CGSolver> Mp_solver;  // For particle solver
    std::shared_ptr<CGSolver> Me_solver;

    std::unique_ptr<mfem::ParGridFunction> CnP;     // Concentration CnP (ParGridFunction)
    std::unique_ptr<mfem::ParGridFunction> CnE;     // Concentration CnE (ParGridFunction)

    // mfem::HypreParMatrix *Mmatp;
    // mfem::HypreParMatrix *Mmate;

    std::shared_ptr<mfem::HypreParMatrix> Mmatp;
    std::shared_ptr<mfem::HypreParMatrix> Mmate;


    mfem::HypreParMatrix *Tmatp;

    // HypreParMatrix Kmatp;
    // HypreParMatrix Kmate;

    std::shared_ptr<mfem::HypreParMatrix> Kmatp;
    std::shared_ptr<mfem::HypreParMatrix> Kmate;

    HypreParVector X1v;
    HypreParVector Fcb;
    HypreParVector Feb;


    std::unique_ptr<mfem::ParGridFunction> Rxc;     
    std::unique_ptr<mfem::ParGridFunction> Rxe;     

    std::unique_ptr<mfem::ParLinearForm> Bc2;
    std::unique_ptr<mfem::ParLinearForm> Be2;

    mfem::ParLinearForm Fct;
    mfem::ParLinearForm Fet;
    // mfem::ParGridFunction PeR;

    // GridFunctionCoefficient matCoef_R;
    // GridFunctionCoefficient cAe;
    GridFunctionCoefficient *cDp;
    GridFunctionCoefficient *cDe;

    std::unique_ptr<mfem::ProductCoefficient> m_nbcCoef;    
    // mfem::Array<int> nbc_w_bdr;       

    double eCrnt;

};

#endif // CONCENTRATIONS_HPP
