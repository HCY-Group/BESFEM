#ifndef REACTION_HPP
#define REACTION_HPP

#include "Concentrations_Base.hpp"

class Reaction {

public:
    Reaction(mfem::ParMesh* pmesh, mfem::ParFiniteElementSpace* fespace, MeshHandler &mh);

    void Initialize(mfem::ParGridFunction &Rx, double initial_value);
    void TimeStep(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn1, mfem::ParGridFunction &Cn2, mfem::ParGridFunction &psx1, mfem::ParGridFunction &psx2, mfem::ParGridFunction &phx1, mfem::ParGridFunction &phx2);
    

private:

    MeshHandler &mesh_handler;

    void SetInitialReaction(mfem::ParGridFunction &Cn, double initial_value);
    void ElectrolyteConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);

    double tc1 =(2*Constants::t_minus-1.0)/(2*Constants::t_minus*(1.0-Constants::t_minus));
	double tc2 = 1.0/(2*Constants::t_minus*(1.0-Constants::t_minus))*Constants::Cst1;

    mfem::ParMesh *pmesh;
    mfem::ParFiniteElementSpace *fespace;

    mfem::ParGridFunction AvP;
    mfem::ParGridFunction AvB;

    int nE;                                         // Number of elements
    int nC;                                         // Number of corners (assuming this is number of corners)
    int nV;                                         // Number of vertices

    mfem::ParGridFunction *Dmp;
    mfem::ParGridFunction *kpl;

    double dffe;


};

#endif // REACTION_HPP