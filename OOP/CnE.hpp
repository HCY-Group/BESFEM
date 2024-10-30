#ifndef CNE_HPP
#define CNE_HPP

#include "Concentrations_Base.hpp"

class CnE : public Concentrations {
public:
    CnE(mfem::ParMesh* pmesh, mfem::ParFiniteElementSpace* fespace, MeshHandler &mh);

protected:
    void Initialize_Differences(mfem::ParGridFunction &psx) override;
};

#endif // CNE_HPP
