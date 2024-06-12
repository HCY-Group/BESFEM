#ifndef MESH_HANDLER_HPP
#define MESH_HANDLER_HPP

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <memory>

class MeshHandler
{
public:
    MeshHandler(const char* mesh_file, const char* dsF_file, int order);
    void InitializeMesh();
    void PrintMeshInfo();

    // Getter functions
    mfem::ParFiniteElementSpace* GetFESpace() const { return fespace.get(); }
    mfem::ParGridFunction* GetPsi() const { return psi.get(); }
    mfem::ParGridFunction* GetPse() const { return pse.get(); }
    double GetGtPsi() const { return gtPsi; }
    mfem::ParMesh* GetPmesh() const { return pmesh.get(); }

private:
    const char* mesh_file;
    const char* dsF_file;
    int order;
    double dh;
    double zeta;
    double eps;
    double rho;
    double Cr;

    mfem::Mesh gmesh;
    std::unique_ptr<mfem::ParMesh> pmesh;
    std::unique_ptr<mfem::FiniteElementSpace> gFespace;
    std::unique_ptr<mfem::ParFiniteElementSpace> fespace;
    std::unique_ptr<mfem::GridFunction> gDsF;
    std::unique_ptr<mfem::ParGridFunction> dsF;
    std::unique_ptr<mfem::ParGridFunction> psi;
    std::unique_ptr<mfem::ParGridFunction> pse;
    std::unique_ptr<mfem::ParGridFunction> AvP;

    void ReadGlobalDistanceFunction();
    void InterpolateDomainParameters();
    void CalculateTotalPsi();
    void CalculateTotalPse();
    //void CalculateTargetCurrent();
    void CalculateTargetCurrent(double tPsi);

    double gtPsi;
    double gtPse;
    double gTrgI;
};

#endif


