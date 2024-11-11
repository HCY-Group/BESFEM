#ifndef POTE_HPP
#define POTE_HPP

#include "Potentials_Base.hpp"

class PotE : public Potentials {

public:
    PotE(mfem::ParMesh* pmesh, mfem::ParFiniteElementSpace* fespace, MeshHandler &mh);
    void Initialize(mfem::ParGridFunction &Cn, double initial_value);

    static mfem::CGSolver *cgPE_solver; // static variable to be used in reaction
    static mfem::CGSolver *GetcgPEsolver() { return cgPE_solver; } // static variable to be used in reaction


private:

    // mfem::CGSolver *cgPE_solver;


};

#endif // POTE_HPP