#ifndef CURRENT_HPP
#define CURRENT_HPP


#include "mfem.hpp"
#include "Mesh_Handler.hpp"

// #include "Concentrations_Base.hpp"

class Current {

public:
    Current(mfem::ParMesh* pmesh, mfem::ParFiniteElementSpace* fespace, MeshHandler &mh);
    
    void Constant(mfem::ParGridFunction &phx, double &global_current);


private:

    MeshHandler &mesh_handler;

    mfem::ParMesh *pmesh;
    mfem::ParFiniteElementSpace *fespace;

    double sgn;
    double dV;



};


#endif // CURRENT_HPP