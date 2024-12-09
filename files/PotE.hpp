#ifndef POTE_HPP
#define POTE_HPP

#include "Potentials_Base.hpp"

extern double BvE;


class PotE : public Potentials {

    std::unique_ptr<mfem::ParBilinearForm> Kl1; // Member variable
    std::unique_ptr<mfem::ParBilinearForm> Kl2; // Member variable



public:
    PotE(mfem::ParMesh* pmesh, mfem::ParFiniteElementSpace* fespace, MeshHandler &mh);
    void Initialize(mfem::ParGridFunction &Cn, double initial_value);
    void TimeStep(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &phx);
    void CalculateGlobalError(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror);


    static mfem::CGSolver *cgPE_solver; // static variable to be used in reaction
    static mfem::CGSolver *GetcgPEsolver() { return cgPE_solver; } // static variable to be used in reaction

    mfem::Array<int> ess_tdof_list_w;
    mfem::ConstantCoefficient dbc_w_Coef;
    mfem::Array<int> dbc_w_bdr;

    double error_E = 1.0;

    // double GetGlobalError() const { return globalerror_E; } // Accessor for external use


private:

    double tc1 =(2*Constants::t_minus-1.0)/(2*Constants::t_minus*(1.0-Constants::t_minus));
	double tc2 = 1.0/(2*Constants::t_minus*(1.0-Constants::t_minus))*Constants::Cst1;

    double dffe;
    
    void ElectrolyteConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);
    // mfem::CGSolver *cgPE_solver;

    mfem::ParGridFunction *Dmp;
    mfem::ParGridFunction *kpl;

    // mfem::ParBilinearForm *Kl1;
    // mfem::ParBilinearForm *Kl2;

    mfem::HypreParMatrix Kdm;
    mfem::HypreParMatrix KmE;

    mfem::ParLinearForm B1t;
    mfem::HypreParVector X1v;
    mfem::HypreParVector B1v;
    mfem::HypreParVector RHSl;


    mfem::HypreParVector *CeVn;
    mfem::HypreParVector *LpCe;

    mfem::HypreSmoother Mpe;

    Array<int> boundary_dofs;

    // mfem::ConstantCoefficient dbc_w_Coef;

    mfem::GridFunctionCoefficient cKp;

    mfem::ParGridFunction *RpE; 
    mfem::ParLinearForm ftPotE; // force term particle electrolyte
    mfem::HypreParVector Flb;
    
    double gtPse;                                   // Total Pse from MeshHandler
    // double globalerror_E; // Member variable to store the global error for PotE


};

#endif // POTE_HPP