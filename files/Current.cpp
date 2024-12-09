#include "Current.hpp"
#include "Mesh_Handler.hpp"
#include "mfem.hpp"
#include "CnE.hpp"
#include "CnP.hpp"
#include "PotE.hpp"
#include "PotP.hpp"
#include "Reaction.hpp"

Current::Current(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh)
    : pmesh(pm), fespace(fe), mesh_handler(mh)

{
    


}

void Current::Constant(mfem::ParGridFunction &phx, double &global_current){

    sgn = copysign(1, gTrgI - global_current);
    dV = Constants::dt * Constants::Vsr * sgn;
    BvP -= dV;
    phx -= dV;

    // std::cout << "BvP: " << BvP << std::endl;

    
}


