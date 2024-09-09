#ifndef POTP_HPP
#define POTP_HPP

#include "mfem.hpp"
#include "MeshHandler.hpp"
#include <memory>

class PotP {

public: 

    PotP(MeshHandler &mesh_handler);

    void Initialize();

    double GetBvP() const { return BvP; }
    const mfem::ParGridFunction& GetphP() const { return phP; } 




private:

    MeshHandler &mesh_handler;

    mfem::ParFiniteElementSpace *fespace;
    mfem::ParGridFunction phP;
    mfem::ParGridFunction kap;
    mfem::ParGridFunction RpP;
    mfem::ParGridFunction pP0;

    mfem::HypreParMatrix KmP;

    double BvP;

};

#endif
