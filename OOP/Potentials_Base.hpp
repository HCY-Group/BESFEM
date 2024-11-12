#ifndef POTENTIALS_HPP
#define POTENTIALS_HPP

// Public: Members can be accessed from anywhere. This is the default access modifier. 
// Protected: Members can be accessed within the class and by classes that inherit from that class. 
// Private: Members can only be accessed within the class that defines them.


#include "mfem.hpp"
#include "Mesh_Handler.hpp"

#include <memory>

class Potentials {
public:
    Potentials(mfem::ParMesh* pmesh, mfem::ParFiniteElementSpace* fespace, MeshHandler &mh);
    virtual ~Potentials() = default;
    MeshHandler &mesh_handler;

    void SetInitialPotentials(mfem::ParGridFunction &ph, double initial_value);
    void SetUpSolver(mfem::CGSolver &solver, double value_1, double value_2);
    void CreateReaction(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, double value);
    void ForceTerm(mfem::ParGridFunction &Rx2, mfem::ParLinearForm &Fxx);

    double BvP = 2.9395;
    double BvE = -1.0;
    double Vcell;

protected:
    
    mfem::ParMesh *pmesh;
    mfem::ParFiniteElementSpace *fespace;


private:

    int nE;                                         // Number of elements
    int nC;                                         // Number of corners (assuming this is number of corners)
    int nV;                                         // Number of vertices

    mfem::ParGridFunction *Rxx;
    mfem::GridFunctionCoefficient *cXx;

};

#endif // POTENTIALS_HPP
