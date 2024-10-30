#ifndef POTENTIALS_HPP
#define POTENTIALS_HPP

#include "mfem.hpp"
#include "Mesh_Handler.hpp"
#include "Constants.hpp"
#include "Reaction.hpp"
#include "Concentrations.hpp"
#include <memory>

class Reaction;

class Potentials {

public:

    Potentials(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh);
    void SetReaction(Reaction *reaction);


    void InitializePotP();
    void InitializePotE();

    mfem::ParGridFunction *phP;
    mfem::ParGridFunction *phE;
    

    mfem::CGSolver *cgPP_solver;
    mfem::CGSolver *cgPE_solver;




private:

    Reaction *reaction;
    
    void CreatePotentials(mfem::ParGridFunction &ph, double initial_value);
    void Solver(mfem::CGSolver &solver, double value_1, double value_2);
    
    MeshHandler &mesh_handler;
    mfem::ParFiniteElementSpace *fespace;
    mfem::ParMesh *pmesh;

    double BvP;
    double BvE;
    double Vcell;

    // int nE;                                         // Number of elements
    // int nC;                                         // Number of corners (assuming this is number of corners)
    // int nV;                                         // Number of vertices

    // mfem::ParGridFunction *phP;
    // mfem::ParGridFunction *phE;

    // mfem::CGSolver *cgPP_solver;
    // mfem::CGSolver *cgPE_solver;

    mfem::HypreParMatrix *KmP; // stiffness matrix PotP
    mfem::HypreParMatrix *KmE; // stiffness matrix PotE
    // mfem::HypreParMatrix *Kdm; // laplace matrix PotE



};



#endif // POTENTIALS_HPP