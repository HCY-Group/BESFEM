#ifndef CNE_HPP
#define CNE_HPP

#include "Concentrations_Base.hpp"

class CnE : public Concentrations {
public:
    CnE(mfem::ParMesh* pmesh, mfem::ParFiniteElementSpace* fespace, MeshHandler &mh);
    void Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx, bool perform_lithiation);

    void TimeStep(mfem::ParGridFunction &Rx);



private:

    mfem::ParGridFunction *PeR;

    std::shared_ptr<mfem::CGSolver> Me_solver;
    std::shared_ptr<mfem::HypreParMatrix> Mmate;
    mfem::HypreSmoother Me_prec;

    mfem::ParGridFunction RxE;

};

#endif // CNE_HPP
