#ifndef REACTION_HPP
#define REACTION_HPP

#include "mfem.hpp"
#include "CnP.hpp"
#include "CnE.hpp"
#include "MeshHandler.hpp"
#include "PotP.hpp"
#include "PotE.hpp"
#include <memory>

class Reaction {
public:
    // Constructor that takes a reference to a MeshHandler object
    Reaction(MeshHandler &mesh_handler, CnE &cne, PotE &pote);

    // Initialize method to set up the initial conditions
    // void Initialize();

    // TimeStep method to perform a single time step
    void TimeStep(double dt);

    // void Save();


private:
    MeshHandler &mesh_handler; // Reference to MeshHandler object
    CnE &cne;
    PotE &pote;

    mfem::ParFiniteElementSpace *fespace; // Pointer to the parallel finite element space
    mfem::ParMesh *pmesh;


    double dffe;
    double tc1;
    double tc2;
    double BvE;

    int nV;
    int vi;


    mfem::ParGridFunction Dmp;
    mfem::ParGridFunction kpl;
    mfem::ParGridFunction CnEGridFunction;
    mfem::ParGridFunction pse; // Grid function for pse
    mfem::ParGridFunction phE;

    mfem::ParLinearForm B1t;

    mfem::HypreParMatrix Kdm;
    mfem::HypreParMatrix Kml;
    mfem::HypreParVector X1v;
    mfem::HypreParVector B1v;
    mfem::HypreParVector CeVn;
    mfem::HypreParVector LpCe;






    // mfem::ParGridFunction pse; // Grid function for psi
    // mfem::ParGridFunction AvP; // Grid function for psi
    // //mfem::ParGridFunction Rxn; // Grid function for psi

    // //std::unique_ptr<mfem::ParGridFunction> Rxn;
    // // mfem::ParGridFunction* Rxn;
    // mfem::ParGridFunction Rxn; // Grid function for Rxn - similar to how psi was in CnP
    // std::unique_ptr<mfem::ParGridFunction> Rxe;
    // // mfem::ParGridFunction Rxe;

    // mfem::ParGridFunction De;

    // double val;
    // double eCrnt;
    // double geCrnt;
    // double infx;
    // double CeC;
    // double gCeC;
    // double CeAvg;

    // int nE;
    // int nC;
    // int nV;
    // int s;

    // double gtPse; // Total Psi value
    // double rho; // Rho value
    // double Cr; // Cr value
    // double L_w;
    // double Ce0;
    // mfem::HypreParMatrix Mmate; // Declare Mmatp as a class member
    // mfem::HypreParMatrix Kmate; // Declare Mmatp as a class member
    // mfem::HypreParMatrix *TmatR; // Declare Mmatp as a class member
    // mfem::HypreParMatrix *TmatL; // Declare Mmatp as a class member
    // mfem::GridFunctionCoefficient matCoef_R; // Declare matCoef_R as a pointer

    // mfem::ParGridFunction CnEGridFunction; // Renamed member variable for CnP grid function
    // mfem::ParGridFunction cDe;

    // mfem::ParGridFunction CeT;

};

#endif // REACTION_HPP