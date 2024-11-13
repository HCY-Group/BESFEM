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
    
    void KMatrix(mfem::ParBilinearForm &K, mfem::GridFunctionCoefficient &gfc, mfem::Array<int> boundary, mfem::ParGridFunction &potential, mfem::ParLinearForm &plf_B, mfem::HypreParMatrix &matrix, mfem::HypreParVector &hpv_X, mfem::HypreParVector &hpv_B);
    void PCG_Solver(mfem::HypreSmoother &smoother, mfem::CGSolver &cg, mfem::HypreParMatrix &KMatrix);

    void CreateReaction(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, double value);
    void ForceTerm(mfem::ParGridFunction &Rx2, mfem::ParLinearForm &Fxx);
    void ForceVector(mfem::ParBilinearForm &K, mfem::Array<int> boundary, mfem::ParGridFunction &potential, mfem::ParLinearForm &plf_B, mfem::HypreParMatrix &matrix, mfem::HypreParVector &hpv_X, mfem::HypreParVector &hpv_B, mfem::ConstantCoefficient &Coef, mfem::Array<int> &bdr);

    double BvP = 2.9395;
    double BvE = -1.0;
    double Vcell;

    int nE;                                         // Number of elements
    int nC;                                         // Number of corners (assuming this is number of corners)
    int nV;                                         // Number of vertices

protected:
    
    mfem::ParMesh *pmesh;
    mfem::ParFiniteElementSpace *fespace;


private:

    // int nE;                                         // Number of elements
    // int nC;                                         // Number of corners (assuming this is number of corners)
    // int nV;                                         // Number of vertices

    mfem::ParGridFunction *Rxx;
    mfem::GridFunctionCoefficient *cXx;

};

#endif // POTENTIALS_HPP
