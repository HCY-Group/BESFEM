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
    mfem::ParGridFunction* GetAvB() const { return AvB.get(); }
    double GetTotalPsi() const { return gtPsi; }
    mfem::ParGridFunction* GetPse() const { return pse.get(); }
    double GetTotalPse() const { return gtPse; }
    double GetLw() const { return L_w; }


private:
    const char* mesh_file;
    const char* dsF_file;
    int order;
    double dh, zeta, eps, rho, Cr;
    double gtPsi, gtPse, gTrgI, L_w, tPsi, val, tPse, trgI;

    int Onm, nV, nE, nC;

    mfem::Mesh gmesh;

    std::unique_ptr<mfem::FiniteElementSpace> gFespace;
    std::unique_ptr<mfem::GridFunction> gDsF;
    std::unique_ptr<mfem::ParMesh> pmesh;
    std::unique_ptr<mfem::ParFiniteElementSpace> fespace;
    std::unique_ptr<mfem::ParGridFunction> dsF, psi, pse, AvP, AvB;

    void ReadGlobalDistanceFunction(const std::unique_ptr<mfem::FiniteElementSpace>& fespace);    
    void CalculateElementVolume(int nE, const std::unique_ptr<mfem::ParMesh>& pmesh, mfem::Vector& EVol);
    void InterpolateDomainParameters(int nV, const std::unique_ptr<mfem::ParFiniteElementSpace>& fespace);
    void CalculateTotalPsi(int nV, int nE, int nC, const mfem::Vector& EVol);
    void CalculateTotalPse(int nV, int nE, int nC, const mfem::Vector& EVol);
    void CalculateTargetCurrent(double tPsi);
};

#endif // MESH_HANDLER_HPP


    // fespace = make_unique<ParFiniteElementSpace>(pmesh.get(), new H1_FECollection(order, pmesh->Dimension()));
    // dsF = make_unique<ParGridFunction>(fespace.get());

    // // Map local distance function from global one
    // Array<HYPRE_BigInt> E_L2G;
    // pmesh->GetGlobalElementIndices(E_L2G);



    // for (int ei = 0; ei < nE; ei++) {
    //     int gei = E_L2G[ei];
    //     gmesh.GetElementVertices(gei, gVTX);
    //     pmesh->GetElementVertices(ei, VTX);
    //     for (int vi = 0; vi < nC; vi++) {
    //         (*dsF)(VTX[vi]) = (*gDsF)(gVTX[vi]);
    //     }
    // }