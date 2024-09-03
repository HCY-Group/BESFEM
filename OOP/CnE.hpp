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

    void Initialize();

    void TimeStep(double dt);

    void Save();

    const mfem::ParGridFunction& GetCnE() const { return CnEGridFunction; } 


private:
    MeshHandler &mesh_handler; 
    CnP &cnp;

    mfem::ParFiniteElementSpace *fespace; 
    mfem::ParGridFunction pse; 
    mfem::ParGridFunction AvP; 

    mfem::ParGridFunction Rxn; // Grid function for Rxn - similar to how psi was in CnP
    std::unique_ptr<mfem::ParGridFunction> Rxe;

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
    int s;

    double gtPse; // Total Psi value
    double rho; // Rho value
    double Cr; // Cr value
    double L_w;
    double Ce0;
    mfem::HypreParMatrix Mmate; 
    mfem::HypreParMatrix Kmate; 
    mfem::HypreParMatrix *TmatR; 
    mfem::HypreParMatrix *TmatL; 
    mfem::GridFunctionCoefficient matCoef_R; 

    mfem::ParGridFunction CnEGridFunction; // Renamed member variable for CnE grid function
    
    
    mfem::ParGridFunction cDe;

    mfem::ParGridFunction CeT;

};

#endif // CNE_HPP
