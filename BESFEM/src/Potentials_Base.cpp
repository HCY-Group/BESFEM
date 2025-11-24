#include "../include/Potentials_Base.hpp"
#include "../include/Initialize_Geometry.hpp"
#include "mfem.hpp"

PotentialBase::PotentialBase(Initialize_Geometry &geo, Domain_Parameters &para)
    : pmesh(geo.parallelMesh.get()), fespace(geo.parfespace), geometry(geo), domain_parameters(para), EVol(para.EVol), px0(fespace.get()),
    TmpF(fespace.get())

{
    nE = geometry.nE; 
    nC = geometry.nC; 
    nV = geometry.nV; 

    X0 = mfem::HypreParVector(fespace.get()); // Initialize the potential vector
    px0 = mfem::ParGridFunction(fespace.get()); // Initialize the potential grid function
    TmpF = mfem::ParGridFunction(fespace.get()); // Temporary grid function for error calculations
}
