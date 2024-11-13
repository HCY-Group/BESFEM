#ifndef POTE_HPP
#define POTE_HPP

#include "Potentials_Base.hpp"

class PotE : public Potentials {

public:
    PotE(mfem::ParMesh* pmesh, mfem::ParFiniteElementSpace* fespace, MeshHandler &mh);
    void Initialize(mfem::ParGridFunction &Cn, double initial_value);
    void TimeStep(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &phx);


    static mfem::CGSolver *cgPE_solver; // static variable to be used in reaction
    static mfem::CGSolver *GetcgPEsolver() { return cgPE_solver; } // static variable to be used in reaction


private:

    double tc1 =(2*Constants::t_minus-1.0)/(2*Constants::t_minus*(1.0-Constants::t_minus));
	double tc2 = 1.0/(2*Constants::t_minus*(1.0-Constants::t_minus))*Constants::Cst1;

    double dffe;
    
    void ElectrolyteConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);
    // mfem::CGSolver *cgPE_solver;

    mfem::ParGridFunction *Dmp;
    mfem::ParGridFunction *kpl;

    mfem::ParBilinearForm *Kl1;
    mfem::ParBilinearForm *Kl2;

    mfem::HypreParMatrix Kdm;
    mfem::HypreParMatrix KmE;

    mfem::ParLinearForm *B1t;
    mfem::HypreParVector *X1v;
    mfem::HypreParVector *B1v;

    mfem::HypreParVector *CeVn;
    mfem::HypreParVector *LpCe;

    mfem::HypreSmoother Mpe;

    Array<int> boundary_dofs;


};

#endif // POTE_HPP