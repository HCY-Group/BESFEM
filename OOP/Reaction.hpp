#ifndef REACTION_HPP
#define REACTION_HPP

#include "Concentrations_Base.hpp"

class Reaction {

public:
    Reaction(mfem::ParMesh* pmesh, mfem::ParFiniteElementSpace* fespace, MeshHandler &mh);

    void Initialize(mfem::ParGridFunction &Rx, double initial_value);

    

private:

    MeshHandler &mesh_handler;

    void SetInitialReaction(mfem::ParGridFunction &Cn, double initial_value);

    mfem::ParMesh *pmesh;
    mfem::ParFiniteElementSpace *fespace;

    mfem::ParGridFunction AvP;
    mfem::ParGridFunction AvB;

    int nE;                                         // Number of elements
    int nC;                                         // Number of corners (assuming this is number of corners)
    int nV;                                         // Number of vertices




};

#endif // REACTION_HPP