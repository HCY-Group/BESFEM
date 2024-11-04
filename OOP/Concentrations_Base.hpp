#ifndef CONCENTRATIONS_HPP
#define CONCENTRATIONS_HPP

#include "mfem.hpp"
#include "Mesh_Handler.hpp"

#include <memory>

class Concentrations {
public:
    Concentrations(mfem::ParMesh* pmesh, mfem::ParFiniteElementSpace* fespace, MeshHandler &mh);
    virtual ~Concentrations() = default;
    MeshHandler &mesh_handler;


    void SetInitialValues(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx, bool perform_lithiation);

    mfem::ParGridFunction psi;                     
    mfem::ParGridFunction pse; 

    // mfem::ParGridFunction *CnP;
    // mfem::ParGridFunction *CnE;

protected:
    mfem::ParMesh *pmesh;
    mfem::ParFiniteElementSpace *fespace;
    std::shared_ptr<mfem::HypreParMatrix> Mmat;
    std::shared_ptr<mfem::CGSolver> solver;
    mfem::HypreSmoother smoother;

    void SetInitialConcentration(mfem::ParGridFunction &Cn, double initial_value);
    void SetUpSolver(mfem::ParGridFunction &psx, std::shared_ptr<mfem::HypreParMatrix> &Mmat, mfem::CGSolver &m_solver, mfem::HypreSmoother &smoother);
    void ImposeNeumannBC(mfem::ParGridFunction &psx, mfem::ParGridFunction &PGF);
    void CreateReaction(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, double value);
    void ForceTerm(mfem::ParGridFunction &gfc, mfem::ParLinearForm &Fxx, mfem::Array<int> boundary, mfem::ProductCoefficient m, bool apply_boundary_conditions);
    void TotalReaction(mfem::ParGridFunction &Rx, double value);
    std::shared_ptr<mfem::GridFunctionCoefficient> Diffusivity(mfem::ParGridFunction &psx, mfem::ParGridFunction &Cn, bool particle_electrolyte );


private:
    // MeshHandler &mesh_handler;

    void LithiationCalculation(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);

    mfem::ParBilinearForm *M = nullptr;
    mfem::ParGridFunction *Ps_gf;
    mfem::GridFunctionCoefficient *cP;

    mfem::Array<int> boundary_dofs;

    mfem::ParGridFunction TmpF;

    int nE;                                         // Number of elements
    int nC;                                         // Number of corners (assuming this is number of corners)
    int nV;                                         // Number of vertices

    const mfem::Vector& EVol;                       // Element volumes from MeshHandler
    double gtPsi;                                   // Total Psi from MeshHandler
    double gtPse;                                   // Total Pse from MeshHandler

    mfem::ParGridFunction *Rxx;
    mfem::GridFunctionCoefficient *cXx;

    double infx;





};

#endif // CONCENTRATIONS_HPP
