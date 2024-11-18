#ifndef POTP_HPP
#define POTP_HPP

#include "Potentials_Base.hpp"

extern double BvP;


class PotP : public Potentials {

    std::unique_ptr<mfem::ParBilinearForm> Kp2; // Member variable


public:
    PotP(mfem::ParMesh* pmesh, mfem::ParFiniteElementSpace* fespace, MeshHandler &mh);
    void Initialize(mfem::ParGridFunction &Cn, double initial_value);
    void TimeStep(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &phx);
    void CalculateGlobalError(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror);


    static mfem::CGSolver *cgPP_solver; // static variable to be used in reaction
    static mfem::CGSolver *GetcgPPsolver() { return cgPP_solver; } // static variable to be used in reaction

    mfem::Array<int> ess_tdof_list_e;
    mfem::ConstantCoefficient dbc_e_Coef;
    mfem::Array<int> dbc_e_bdr;

    double error_P = 1.0;

    // double GetGlobalError() const { return globalerror_P; } // Accessor for external use



private:

    void ParticleConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);


    // mfem::CGSolver *cgPP_solver;

    mfem::ParGridFunction *RpP; 
    mfem::ParLinearForm ftPotP; // force term particle potential
    mfem::HypreParVector Fpb;

    mfem::ParGridFunction *kap;
    // std::unique_ptr<mfem::ParBilinearForm> Kp2; // Change to unique_ptr
    mfem::ParLinearForm B1t;
    mfem::HypreParVector X1v;
    mfem::HypreParVector B1v;
    mfem::HypreSmoother Mpp;
    mfem::HypreParMatrix KmP;

    double gtPsi;                                   // Total Psi from MeshHandler
    // double error_P;
    // double globalerror_P;

    // double globalerror_P; // Member variable to store the global error for PotP



};

#endif // POTP_HPP