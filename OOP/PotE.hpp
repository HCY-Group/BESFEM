#ifndef POTE_HPP
#define POTE_HPP

#include "Potentials_Base.hpp"

class PotE : public Potentials {

public:
    PotE(mfem::ParMesh* pmesh, mfem::ParFiniteElementSpace* fespace, MeshHandler &mh);
    void Initialize(mfem::ParGridFunction &Cn, double initial_value);


private:

    mfem::CGSolver *cgPE_solver;


};

#endif // POTE_HPP