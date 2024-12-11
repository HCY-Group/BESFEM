#ifndef CNE_HPP
#define CNE_HPP

#include "Concentrations_Base.hpp"

class CnE : public Concentrations {
public:
    CnE(mfem::ParMesh* pmesh, mfem::ParFiniteElementSpace* fespace, MeshHandler &mh);
    void Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx);

    void TimeStep(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);
    mfem::ParGridFunction *RxE;
    mfem::ParLinearForm ftE;
    mfem::Array<int> nbc_w_bdr;
    std::unique_ptr<mfem::ProductCoefficient> m_nbcCoef; 
    mfem::Array<int> boundary_dofs;

    static mfem::HypreParVector *CeVn; // static variable to be used in reaction
    static mfem::HypreParVector* GetCeVn() { return CeVn; } // static variable to be used in reaction

private:


    mfem::ParGridFunction *PeR;

    std::shared_ptr<mfem::CGSolver> Me_solver;
    std::shared_ptr<mfem::HypreParMatrix> Mmate;
    mfem::HypreSmoother Me_prec;

    // mfem::ParGridFunction RxE;

    double eCrnt;
    // double infx;
    // mfem::ParLinearForm ftE;

    // mfem::Array<int> nbc_w_bdr;
    // std::unique_ptr<mfem::ProductCoefficient> m_nbcCoef;    

    std::shared_ptr<mfem::HypreParMatrix> Kmate;

    mfem::HypreParVector Feb;
    // mfem::Array<int> boundary_dofs;
    mfem::HypreParVector X1v;

    mfem::HypreParVector *CeV0;
    mfem::HypreParVector *RHCe;
    mfem::HypreParMatrix *TmatR;
    mfem::HypreParMatrix *TmatL;



};

#endif // CNE_HPP
