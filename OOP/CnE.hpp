#ifndef CNE_HPP
#define CNE_HPP

#include "Concentrations_Base.hpp"

class CnE : public Concentrations {
public:
    CnE(mfem::ParMesh* pmesh, mfem::ParFiniteElementSpace* fespace, MeshHandler &mh);
    void Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx, bool perform_lithiation);


protected:
    // void Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx, std::shared_ptr<mfem::HypreParMatrix> &Mmat, mfem::CGSolver &m_solver, mfem::HypreSmoother &smoother, bool perform_lithiation);



private:

    mfem::ParGridFunction *PeR;

    std::shared_ptr<mfem::CGSolver> Me_solver;
    std::shared_ptr<mfem::HypreParMatrix> Mmate;
    mfem::HypreSmoother Me_prec;
};

#endif // CNE_HPP
