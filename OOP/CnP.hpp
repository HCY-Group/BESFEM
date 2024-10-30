#ifndef CNP_HPP
#define CNP_HPP

#include "Concentrations_Base.hpp"

class CnP : public Concentrations {

public:
    CnP(mfem::ParMesh* pmesh, mfem::ParFiniteElementSpace* fespace, MeshHandler &mh);
    void Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx, bool perform_lithiation);

    
protected:
    // void Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx, std::shared_ptr<mfem::HypreParMatrix> &Mmat, mfem::CGSolver &m_solver, mfem::HypreSmoother &smoother, bool perform_lithiation);


private:

    mfem::HypreParVector PsVc;

    std::shared_ptr<mfem::CGSolver> Mp_solver;
    std::shared_ptr<mfem::HypreParMatrix> Mmatp;
    mfem::HypreSmoother Mp_prec;


};

#endif // CNP_HPP
