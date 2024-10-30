#ifndef CNP_HPP
#define CNP_HPP

#include "Concentrations_Base.hpp"

class CnP : public Concentrations {
public:
    CnP(mfem::ParMesh* pmesh, mfem::ParFiniteElementSpace* fespace, MeshHandler &mh);
    
protected:
    void Initialize_Differences(mfem::ParGridFunction &psx) override;


private:

    mfem::HypreParVector PsVc;


};

#endif // CNP_HPP
