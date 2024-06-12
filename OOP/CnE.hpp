#ifndef CNE_HPP
#define CNE_HPP

#include "mfem.hpp"
#include <memory>

void InitializeCnE(mfem::ParFiniteElementSpace &fespace,
                mfem::ParGridFunction &pse,
                mfem::ParMesh &pmesh);

#endif
