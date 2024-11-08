#ifndef POTP_HPP
#define POTP_HPP

#include "Potentials_Base.hpp"

class PotP : public Potentials {

public:
    PotP(mfem::ParMesh* pmesh, mfem::ParFiniteElementSpace* fespace, MeshHandler &mh);
    void Initialize(mfem::ParGridFunction &Cn, double initial_value);


private:

    mfem::CGSolver *cgPP_solver;


};

#endif // POTP_HPP