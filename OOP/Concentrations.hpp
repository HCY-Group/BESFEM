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
    Concentrations(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh);

    // Initialization method
    void InitializeCnP();
    void InitializeCnE();

    void TimeStepCnP();
    void TimeStepCnE();

    void SaveCnP();
    void SaveCnE();


    // void SetupBoundaryConditions(std::shared_ptr<ParFiniteElementSpace> fespace);

    // void TestFESpace(ParFiniteElementSpace *fespace);

    Array<int> nbc_w_bdr;

    Reaction reaction;

    // mfem::ParFiniteElementSpace* GetFESpace() const { return fespace.get(); }

    // std::shared_ptr<ParFiniteElementSpace> fespace;
    // void TestFESpace(std::shared_ptr<ParFiniteElementSpace> fespace);


private:

    MeshHandler &mesh_handler;
    // Reaction reaction;

    // std::shared_ptr<mfem::ParMesh> pmesh;
    // std::shared_ptr<mfem::ParFiniteElementSpace> fespace;

    mfem::ParFiniteElementSpace *fespace;
    mfem::ParMesh *pmesh;

    // Reaction reaction;

    // Internal methods for different operations
    void CreateCnE(mfem::ParGridFunction &Cn, double initial_value);
    void Lithiation(mfem::ParGridFunction &Cn, double initial_value);
    void LithiationCalculation(mfem::ParGridFunction &Cn);
    void Solver(mfem::ParGridFunction &psx, std::shared_ptr<mfem::HypreParMatrix> &Mmat, mfem::CGSolver &m_solver, mfem::HypreSmoother &smoother);
    void ImposeNeumannBC(mfem::ParGridFunction &psx, mfem::ParGridFunction &PGF);
    // std::shared_ptr<mfem::GridFunctionCoefficient> ImposeNeumannBC(mfem::ParGridFunction &psx);

    void SetupRx(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, double value, GridFunctionCoefficient cAx);
    void ForceTerm(mfem::ParGridFunction &gfc, mfem::ParLinearForm &Fxx, Array<int> boundary, ProductCoefficient m, bool apply_boundary_conditions);
    void TotalReaction(mfem::ParGridFunction &Rx, double xCrnt);
    void EnsureValidBoundaryAndFESpace();
    void DebugBoundaryArray(const Array<int> &boundary);
    std::shared_ptr<mfem::GridFunctionCoefficient> Diffusivity(mfem::ParGridFunction &psx, mfem::ParGridFunction &Cn, bool particle_electrolyte );
    void K_Matrix(Array<int> boundary, mfem::ParGridFunction &Cn, ParLinearForm &Fxx, std::shared_ptr<HypreParMatrix> &Kmatx, HypreParVector &X1v, HypreParVector &Fxb, GridFunctionCoefficient *cDx);
    void SaltConservation(mfem::ParGridFunction &Cn);


    // Member variables
    // mfem::ParFiniteElementSpace* fespace;           // Finite element space

    // std::shared_ptr<mfem::ParFiniteElementSpace> fespace;

    mfem::ParGridFunction psi;                     // Psi grid function from MeshHandler
    mfem::ParGridFunction pse;                     // Pse grid function from MeshHandler
    mfem::ParGridFunction TmpF;
    // mfem::ParMesh* pmesh;

    const mfem::Vector& EVol;                       // Element volumes from MeshHandler
    double gtPsi;                                   // Total Psi from MeshHandler
    double gtPse;                                   // Total Pse from MeshHandler
    double infx;

    int nE;                                         // Number of elements
    int nC;                                         // Number of corners (assuming this is number of corners)
    int nV;                                         // Number of vertices

    Array<int> boundary_dofs;

    // std::shared_ptr<CGSolver> Mp_solver;  // For particle solver
    // std::shared_ptr<CGSolver> Me_solver;
    // std::shared_ptr<CGSolver> solver;

    mfem::HypreSmoother Mp_prec;
    mfem::HypreSmoother Me_prec;

    mfem::CGSolver *m_solver = nullptr;
    mfem::CGSolver *solver;
    // mfem::CGSolver *Mp_solver;
    std::shared_ptr<mfem::CGSolver> Mp_solver;
    std::shared_ptr<mfem::CGSolver> Me_solver;

    mfem::ParGridFunction *CnP;
    mfem::ParGridFunction *CnE;

    // mfem::HypreParMatrix *Mmatp;
    // mfem::HypreParMatrix *Mmate;

    mfem::HypreParMatrix *Mmat;

    std::shared_ptr<mfem::HypreParMatrix> Mmatp;
    std::shared_ptr<mfem::HypreParMatrix> Mmate;    
    
    mfem::HypreParMatrix *Tmatp;
    mfem::HypreParMatrix *TmatR;
    mfem::HypreParMatrix *TmatL;

    // HypreParMatrix Kmatp;
    // HypreParMatrix Kmate;

    std::shared_ptr<mfem::HypreParMatrix> Kmatp;
    std::shared_ptr<mfem::HypreParMatrix> Kmate;

    HypreParVector X1v;
    HypreParVector Fcb;
    HypreParVector Feb;

    // mfem::HypreParVector *Fv_ptr;
    // mfem::HypreParVector F_hpv;


    // std::unique_ptr<mfem::ParGridFunction> Rxc;     
    // std::unique_ptr<mfem::ParGridFunction> Rxe;   
    
    mfem::ParGridFunction *Rxn;
    mfem::ParGridFunction *Rxc; 
    mfem::ParGridFunction *Rxe; 

    mfem::ParGridFunction *PGF;
    mfem::ParGridFunction *PeR;

    mfem::GridFunctionCoefficient *cXx;
    mfem::GridFunctionCoefficient *matCoef_R;
    std::shared_ptr<mfem::GridFunctionCoefficient> coef;
    
    // mfem::ParLinearForm *Bx2; 

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

    double Ce0 = 0.001;
    double CeC = 0.0;
    double gCeC = 0.0;
    double CeAvg = 0.0;

    mfem::ParBilinearForm *M = nullptr;
    mfem::GridFunctionCoefficient *cP;
    mfem::ParGridFunction *Ps_gf;
    mfem::ParGridFunction *Rxx;

    mfem::HypreParVector PsVc;

};

#endif // CONCENTRATIONS_HPP
