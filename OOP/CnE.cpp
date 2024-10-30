#include "CnE.hpp"
#include "mfem.hpp"


CnE::CnE(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh)
    : Concentrations(pm, fe, mh) {}

void CnE::Initialize_Differences(mfem::ParGridFunction &psx) {

    // Impose Neumann boundary condition for CnE
    mfem::ParGridFunction PGF(fespace);
    ImposeNeumannBC(psx, PGF);

    // PGF.Print(std::cout);

}
