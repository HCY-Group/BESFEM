#ifndef CNP_HPP
#define CNP_HPP

#include "mfem.hpp"
#include <memory>

void InitializeCnP(mfem::ParFiniteElementSpace &fespace, 
                       mfem::ParGridFunction &psi, 
                       double gtPsi, mfem::ParMesh &pmesh);

#endif
