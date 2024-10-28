// Reaction.hpp
#ifndef REACTION_HPP
#define REACTION_HPP

#include "mfem.hpp"
#include "Mesh_Handler.hpp"
#include "Constants.hpp"
#include "Concentrations.hpp"

#include <memory>

using namespace mfem;

class Potentials;
class Concentrations; // Forward declaration

class Reaction {
public:
    Reaction(mfem::ParFiniteElementSpace *fe, MeshHandler &mh, Concentrations &con); // Constructor
    void SetPotentials(Potentials *pot);

    void Initialize(); // Method to initialize Rxn
    void TimeStep();

    // std::unique_ptr<mfem::ParGridFunction> Rxn;     // Rxn (ParGridFunction)

    // mfem::ParGridFunction* GetRxn() const { return Rxn.get(); }

    mfem::ParGridFunction *Rxn;

private:
    
    MeshHandler &mesh_handler;
    Concentrations &concentrations; // Reference to Concentrations
    Potentials *potentials;
    // ParGridFunction Rxn; // ParGridFunction for Rxn

    void CreateRx(mfem::ParGridFunction &Rx, double initial_value);
    void SetAvP(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Av, double value);
    void ElectrolyteConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);
    void KMatrix(mfem::ParBilinearForm &K, mfem::GridFunctionCoefficient &gfc, mfem::Array<int> boundary, mfem::ParGridFunction &potential, mfem::ParLinearForm &plf_B, mfem::HypreParMatrix &matrix, mfem::HypreParVector &hpv_X, mfem::HypreParVector &hpv_B);
    void PCG_Solver(mfem::HypreSmoother &smoother, mfem::CGSolver &cg, mfem::HypreParMatrix &KMatrix);
    void ParticleConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);
    void ExchangeCurrentDensity(mfem::ParGridFunction &Av_pgf, mfem::ParGridFunction &Cn);
    void SetZero(mfem::ParGridFunction &pgf);

    int nV;                                         // Number of vertices
    int val;
    int inlp = 0;

    Array<int> boundary_dofs;

    double dffe; 
    double tc1 =(2*Constants::t_minus-1.0)/(2*Constants::t_minus*(1.0-Constants::t_minus));
	double tc2 = 1.0/(2*Constants::t_minus*(1.0-Constants::t_minus))*Constants::Cst1;

    double gErrP = 1.0;
    double gErrE = 1.0;

    mfem::ParGridFunction *CnP;
    mfem::ParGridFunction *CnE;
    mfem::ParFiniteElementSpace *fespace;

    mfem::ParGridFunction psi;                     
    mfem::ParGridFunction pse; 

    mfem::ParGridFunction *phE;                     
    mfem::ParGridFunction *phP;                    

    mfem::ParGridFunction AvP;
    mfem::ParGridFunction *AvP_PGF;
    mfem::ParGridFunction AvB;
    mfem::ParGridFunction *AvB_PGF;

    mfem::ParGridFunction *Dmp;
    mfem::ParGridFunction *kpl;
    mfem::ParGridFunction *kap;


    mfem::ParBilinearForm *Kl1;
    mfem::ParLinearForm *B1t;
    mfem::HypreParVector *X1v;
    mfem::HypreParVector *B1v;

    mfem::HypreParMatrix Kdm;
    
    mfem::HypreParVector CeVn;
    mfem::HypreParVector *LpCe;

    mfem::ParBilinearForm *Kl2;
    mfem::HypreParMatrix KmE;

    mfem::HypreSmoother Mpp;
    mfem::HypreSmoother Mpe;


    mfem::ParBilinearForm *Kp2;
    mfem::HypreParMatrix KmP;


    mfem::ParGridFunction *i0C;
    mfem::ParGridFunction *OCV;
    mfem::ParGridFunction *Kfw;
    mfem::ParGridFunction *Kbw;







};

#endif // REACTION_HPP
