#ifndef POTP_HPP
#define POTP_HPP

#include "Potentials_Base.hpp"

class PotP : public Potentials {

public:
    PotP(mfem::ParMesh* pmesh, mfem::ParFiniteElementSpace* fespace, MeshHandler &mh);
    void Initialize(mfem::ParGridFunction &Cn, double initial_value);
    void CalculateGlobalError(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx);


    static mfem::CGSolver *cgPP_solver; // static variable to be used in reaction
    static mfem::CGSolver *GetcgPPsolver() { return cgPP_solver; } // static variable to be used in reaction



private:

    // mfem::CGSolver *cgPP_solver;

    mfem::ParGridFunction *RpP; 
    mfem::ParLinearForm ftPotP; // force term particle potential
    mfem::HypreParVector Fpb;


};

#endif // POTP_HPP