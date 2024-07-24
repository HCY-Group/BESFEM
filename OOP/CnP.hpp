#ifndef CNP_HPP
#define CNP_HPP

#include "mfem.hpp"
#include "MeshHandler.hpp"
#include <memory>

class CnP {
public:
    // Constructor that takes a reference to a MeshHandler object
    CnP(MeshHandler &mesh_handler);

    // Initialize method to set up the initial conditions
    void Initialize();

    // TimeStep method to perform a single time step
    void TimeStep(double dt);

    void Save();

    mfem::ParGridFunction* GetRxn() const { return Rxn.get(); }



private:
    MeshHandler &mesh_handler; // Reference to MeshHandler object

    mfem::ParFiniteElementSpace *fespace; // Pointer to the parallel finite element space
    mfem::ParGridFunction psi; // Grid function for psi
    mfem::ParGridFunction AvP; // Grid function for psi
    
   //mfem::ParGridFunction Rxn; // Grid function for psi

    std::unique_ptr<mfem::ParGridFunction> Rxn;



    double gtPsi; // Total Psi value
    double rho; // Rho value
    double Cr; // Cr value
    mfem::HypreParMatrix Mmatp; // Declare Mmatp as a class member
    mfem::HypreParMatrix Kmatp; // Declare Mmatp as a class member
    mfem::HypreParMatrix *Tmatp; // Declare Mmatp as a class member

    mfem::ParGridFunction CnPGridFunction; // Renamed member variable for CnP grid function

};

#endif // CNP_HPP
