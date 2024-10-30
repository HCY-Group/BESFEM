#include "CnP.hpp"
#include "mfem.hpp"


CnP::CnP(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh)
    : Concentrations(pm, fe, mh) 
    
    {

    PsVc = mfem::HypreParVector(fespace);


    }

void CnP::Initialize_Differences(mfem::ParGridFunction &psx) {

    psx.GetTrueDofs(PsVc);

}
