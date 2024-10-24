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


    int nV;                                         // Number of vertices

    Array<int> boundary_dofs;

    double dffe; 
    double tc1 =(2*Constants::t_minus-1.0)/(2*Constants::t_minus*(1.0-Constants::t_minus));
	double tc2 = 1.0/(2*Constants::t_minus*(1.0-Constants::t_minus))*Constants::Cst1;

    mfem::ParGridFunction *CnP;
    mfem::ParGridFunction *CnE;
    mfem::ParFiniteElementSpace *fespace;

    mfem::ParGridFunction psi;                     
    mfem::ParGridFunction pse; 

    mfem::ParGridFunction *phE;                     
    mfem::ParGridFunction *phP;                    

    mfem::ParGridFunction AvP;
    mfem::ParGridFunction *AvP_PGF;

    mfem::ParGridFunction *Dmp;
    mfem::ParGridFunction *kpl;

    mfem::ParBilinearForm *Kl1;
    mfem::ParLinearForm *B1t;
    mfem::HypreParVector *X1v;
    mfem::HypreParVector *B1v;

    mfem::HypreParMatrix Kdm;
    
    mfem::HypreParVector CeVn;
    mfem::HypreParVector *LpCe;








};

#endif // REACTION_HPP
