#ifndef REACTION_HPP
#define REACTION_HPP

#include "Concentrations_Base.hpp"

class Reaction {

public:
    Reaction(mfem::ParMesh* pmesh, mfem::ParFiniteElementSpace* fespace, MeshHandler &mh);

    void Initialize(mfem::ParGridFunction &Rx, double initial_value);
    void ExchangeCurrentDensity(mfem::ParGridFunction &Cn);
    void ButlerVolmer(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn1, mfem::ParGridFunction &Cn2, mfem::ParGridFunction &phx1, mfem::ParGridFunction &phx2);
    void TotalReactionCurrent(mfem::ParGridFunction &Rx, double &global_current);

    double global_current;


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

    double local_current;
    // double global_current;

    mfem::ParGridFunction *i0C;
    mfem::ParGridFunction *OCV;
    mfem::ParGridFunction *Kfw;
    mfem::ParGridFunction *Kbw;
    mfem::ParGridFunction *dPHE;

    const mfem::Vector& EVol;                       // Element volumes from MeshHandler






    // double dffe;


};

#endif // REACTION_HPP