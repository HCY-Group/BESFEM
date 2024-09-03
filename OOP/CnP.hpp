#ifndef CNP_HPP
#define CNP_HPP

#include "mfem.hpp"
#include "MeshHandler.hpp"
#include <memory>

class CnP {
public:
    // Constructor that takes a reference to a MeshHandler object
    CnP(MeshHandler &mesh_handler);

    void Initialize();

    void TimeStep(double dt);

    void Save();

    mfem::ParGridFunction* GetRxn() const { return Rxn.get(); } // similar to how psi came from meshhandler
    // .get() needed since Rxn is a unique ptr (helps memory management)



private:
    MeshHandler &mesh_handler; 

    mfem::ParFiniteElementSpace *fespace; 
    mfem::ParGridFunction psi; 
    mfem::ParGridFunction AvP; 
    
    std::unique_ptr<mfem::ParGridFunction> Rxn; // similar to how psi was done

    int nE;
    int nC;
    int nV;
    int nDof;
    
    double Cp0;
    double Xfr;
    double lSum;
    double val;
    double gSum;

    double gtPsi; // Total Psi value
    double rho; // Rho value
    double Cr; // Cr value
    mfem::HypreParMatrix Mmatp; 
    mfem::HypreParMatrix Kmatp; 
    mfem::HypreParMatrix *Tmatp; 

    mfem::ParGridFunction CnPGridFunction; // Renamed member variable for CnP grid function

};

#endif // CNP_HPP
