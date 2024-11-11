#ifndef REACTION_HPP
#define REACTION_HPP

#include "Concentrations_Base.hpp"

class Reaction {

public:
    Reaction(mfem::ParMesh* pmesh, mfem::ParFiniteElementSpace* fespace, MeshHandler &mh);

    void Initialize(mfem::ParGridFunction &Rx, double initial_value);
    void TimeStep(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn1, mfem::ParGridFunction &Cn2, mfem::ParGridFunction &psx1, mfem::ParGridFunction &psx2, mfem::ParGridFunction &phx1, mfem::ParGridFunction &phx2);
    void ButlerVolmer(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn1, mfem::ParGridFunction &Cn2, mfem::ParGridFunction &phx1, mfem::ParGridFunction &phx2);


private:

    MeshHandler &mesh_handler;

    void SetInitialReaction(mfem::ParGridFunction &Cn, double initial_value);
    void ElectrolyteConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);
    void ParticleConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);
    void KMatrix(mfem::ParBilinearForm &K, mfem::GridFunctionCoefficient &gfc, mfem::Array<int> boundary, mfem::ParGridFunction &potential, mfem::ParLinearForm &plf_B, mfem::HypreParMatrix &matrix, mfem::HypreParVector &hpv_X, mfem::HypreParVector &hpv_B);
    void PCG_Solver(mfem::HypreSmoother &smoother, mfem::CGSolver &cg, mfem::HypreParMatrix &KMatrix);
    void ExchangeCurrentDensity(mfem::ParGridFunction &Av_pgf, mfem::ParGridFunction &Cn);
    // void ButlerVolmer(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn1, mfem::ParGridFunction &Cn2, mfem::ParGridFunction &phx1, mfem::ParGridFunction &phx2);

    double tc1 =(2*Constants::t_minus-1.0)/(2*Constants::t_minus*(1.0-Constants::t_minus));
	double tc2 = 1.0/(2*Constants::t_minus*(1.0-Constants::t_minus))*Constants::Cst1;

    Array<int> boundary_dofs;


    mfem::ParMesh *pmesh;
    mfem::ParFiniteElementSpace *fespace;

    mfem::ParGridFunction AvP;
    mfem::ParGridFunction AvB;

    int nE;                                         // Number of elements
    int nC;                                         // Number of corners (assuming this is number of corners)
    int nV;                                         // Number of vertices

    mfem::ParGridFunction *Dmp;
    mfem::ParGridFunction *kpl;
    mfem::ParGridFunction *kap;


    mfem::ParBilinearForm *Kl1;
    mfem::ParBilinearForm *Kl2;
    mfem::ParBilinearForm *Kp2;

    mfem::ParLinearForm *B1t;
    mfem::HypreParVector *X1v;
    mfem::HypreParVector *B1v;

    mfem::HypreParMatrix Kdm;
    mfem::HypreParMatrix KmE;
    mfem::HypreParMatrix KmP;

    mfem::HypreParVector CeVn;
    mfem::HypreParVector *LpCe;

    mfem::HypreSmoother Mpp;
    mfem::HypreSmoother Mpe;

    mfem::CGSolver cgPE_solver;
    mfem::CGSolver cgPP_solver;

    mfem::ParGridFunction *i0C;
    mfem::ParGridFunction *OCV;
    mfem::ParGridFunction *Kfw;
    mfem::ParGridFunction *Kbw;
    mfem::ParGridFunction *dPHE;






    double dffe;


};

#endif // REACTION_HPP