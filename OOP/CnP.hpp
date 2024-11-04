#ifndef CNP_HPP
#define CNP_HPP

#include "Concentrations_Base.hpp"

class CnP : public Concentrations {

public:
    CnP(mfem::ParMesh* pmesh, mfem::ParFiniteElementSpace* fespace, MeshHandler &mh);
    void Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx, bool perform_lithiation);

    void TimeStep(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);

    mfem::ParGridFunction *RxP;

    

private:

    mfem::HypreParVector PsVc;

    std::shared_ptr<mfem::CGSolver> Mp_solver;
    std::shared_ptr<mfem::HypreParMatrix> Mmatp;
    mfem::HypreSmoother Mp_prec;

    mfem::ParLinearForm ftP;


};

#endif // CNP_HPP
