#ifndef CONCENTRATIONS_HPP
#define CONCENTRATIONS_HPP

#include "mfem.hpp"
#include "Mesh_Handler.hpp"
#include "Constants.hpp"
#include <memory>

class Concentrations {
public:
    // Constructor that takes a MeshHandler reference
    Concentrations(MeshHandler &mesh_handler);

    // Initialization method
    void InitializeCnP();

private:
    // Internal methods for different operations
    void Lithiation(mfem::ParGridFunction &Cn, double initial_value);
    void SBM_Matrix(mfem::ParGridFunction &psx, HypreParMatrix &Mmat);
    void Solver(HypreParMatrix &Mmat);
    void SetupBoundaryConditions();

    // Member variables
    mfem::ParFiniteElementSpace* fespace;           // Finite element space
    mfem::ParGridFunction psi;                     // Psi grid function from MeshHandler
    mfem::ParGridFunction pse;                     // Pse grid function from MeshHandler
    mfem::ParGridFunction TmpF;

    const mfem::Vector& EVol;                       // Element volumes from MeshHandler
    double gtPsi;                                   // Total Psi from MeshHandler

    int nE;                                         // Number of elements
    int nC;                                         // Number of corners (assuming this is number of corners)
    int nV;                                         // Number of vertices

    Array<int> boundary_dofs;

    std::unique_ptr<mfem::ParGridFunction> CnP;     // Concentration CnP (ParGridFunction)
    std::unique_ptr<mfem::ParGridFunction> CnE;     // Concentration CnE (ParGridFunction)

    HypreParMatrix Mmatp;
    HypreParMatrix Mmate;
};

#endif // CONCENTRATIONS_HPP
