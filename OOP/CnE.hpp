#ifndef CNE_HPP
#define CNE_HPP

#include "mfem.hpp"
#include "CnP.hpp"
#include "MeshHandler.hpp"
#include <memory>

class CnE {
public:
    // Constructor that takes a reference to a MeshHandler object
    CnE(MeshHandler &mesh_handler, CnP &cnp);

    // Initialize method to set up the initial conditions
    void Initialize();

    // TimeStep method to perform a single time step
    void TimeStep(double dt);

    void Save();


private:
    MeshHandler &mesh_handler; // Reference to MeshHandler object
    CnP &cnp;

    mfem::ParFiniteElementSpace *fespace; // Pointer to the parallel finite element space
    mfem::ParGridFunction pse; // Grid function for psi
    mfem::ParGridFunction AvP; // Grid function for psi
    //mfem::ParGridFunction Rxn; // Grid function for psi

    //std::unique_ptr<mfem::ParGridFunction> Rxn;
    // mfem::ParGridFunction* Rxn;
    mfem::ParGridFunction Rxn; // Grid function for Rxn - similar to how psi was in CnP
    std::unique_ptr<mfem::ParGridFunction> Rxe;
    // mfem::ParGridFunction Rxe;

    mfem::ParGridFunction De;

    double val;
    double eCrnt;
    double geCrnt;
    double infx;
    double CeC;
    double gCeC;
    double CeAvg;

    int nE;
    int nC;
    int nV;

    double gtPse; // Total Psi value
    double rho; // Rho value
    double Cr; // Cr value
    double L_w;
    double Ce0;
    mfem::HypreParMatrix Mmate; // Declare Mmatp as a class member
    mfem::HypreParMatrix Kmate; // Declare Mmatp as a class member
    mfem::HypreParMatrix *TmatR; // Declare Mmatp as a class member
    mfem::HypreParMatrix *TmatL; // Declare Mmatp as a class member
    mfem::GridFunctionCoefficient matCoef_R; // Declare matCoef_R as a pointer

    mfem::ParGridFunction CnEGridFunction; // Renamed member variable for CnP grid function
    mfem::ParGridFunction cDe;

    mfem::ParGridFunction CeT;

};

#endif // CNE_HPP
