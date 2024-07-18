#ifndef MESH_HANDLER_HPP
#define MESH_HANDLER_HPP

#include "mfem.hpp"

class MeshHandler {
public:
    MeshHandler();
    void LoadMesh();
    void InitializeMesh();
    void PrintMeshInfo();
    void Save();
    
    mfem::ParMesh* GetParMesh() const { return pmesh.get(); }
    mfem::ParFiniteElementSpace* GetFESpace() const { return fespace.get(); }
    mfem::ParGridFunction* GetPsi() const { return psi.get(); }
    mfem::ParGridFunction* GetAvP() const { return AvP.get(); }
    double GetTotalPsi() const { return gtPsi; }
    mfem::ParGridFunction* GetPse() const { return pse.get(); }
    double GetTotalPse() const { return gtPse; }
    double GetLw() const { return L_w; }


private:
    const char* mesh_file;
    const char* dsF_file;
    int order;
    double dh, zeta, eps, rho, Cr;
    double gtPsi, gtPse, gTrgI, L_w;

    mfem::Mesh gmesh;
    std::unique_ptr<mfem::FiniteElementSpace> gFespace;
    std::unique_ptr<mfem::GridFunction> gDsF;
    std::unique_ptr<mfem::ParMesh> pmesh;
    std::unique_ptr<mfem::ParFiniteElementSpace> fespace;
    std::unique_ptr<mfem::ParGridFunction> dsF, psi, pse, AvP;

    void ReadGlobalDistanceFunction();
    void InterpolateDomainParameters();
    void CalculateTotalPsi();
    void CalculateTotalPse();
    void CalculateTargetCurrent(double tPsi);
};

#endif // MESH_HANDLER_HPP
