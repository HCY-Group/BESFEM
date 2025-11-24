#include "../include/Concentrations_Base.hpp"
#include "../inputs/Constants.hpp"
#include "../include/Initialize_Geometry.hpp"
#include "mfem.hpp"

ConcentrationBase::ConcentrationBase(Initialize_Geometry &geo, Domain_Parameters &para)
    : pmesh(geo.parallelMesh.get()), fespace(geo.parfespace), geometry(geo), domain_parameters(para), EVol(para.EVol), gtPsi(para.gtPsi), gtPse(para.gtPse),
    nE(geo.nE), nC(geo.nC), nV(geo.nV), VtxVal(geo.nC), EAvg(geo.nE), gmesh(geo.globalMesh.get())


{
    // Allocate once for reuse
    VtxVal.SetSize(nC); // Set size for nodal values
    EAvg.SetSize(nE);   // Set size for average element contributions
    TmpF = std::make_unique<mfem::ParGridFunction>(fespace.get());

}
