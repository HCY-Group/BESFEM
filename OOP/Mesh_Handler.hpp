#ifndef MESH_HANDLER_HPP
#define MESH_HANDLER_HPP

#include "Constants.hpp"
#include "mfem.hpp"
#include <memory>

using namespace mfem;

class MeshHandler {
public:
    // Constructor
    MeshHandler();

    // Public member functions
    void LoadMesh();
    void Save();
    void SetupBoundaryConditions(mfem::ParMesh *pmesh, mfem::ParFiniteElementSpace *fespace);
    // void TestFESpace();

    // Functions to Get Values to Use Elsewhere
    // mfem::ParFiniteElementSpace* GetFESpace() const { return fespace.get(); }

    std::shared_ptr<mfem::ParFiniteElementSpace> GetFESpace();


    mfem::ParMesh* GetPmesh() const { return pmesh.get(); }
    mfem::ParGridFunction* GetPsi() const { return psi.get(); }
    mfem::ParGridFunction* GetPse() const { return pse.get(); }
    mfem::ParGridFunction* GetAvP() const { return AvP.get(); }
    mfem::ParGridFunction* GetAvB() const { return AvB.get(); }


    mfem::ParMesh GetMesh();
    

    // const mfem::Vector& GetElementVolume() const {return EVol; }
    // const mfem::Vector& EVol = mesh_handler.GetElementVolume();
    const mfem::Vector& GetElementVolume() const;

    std::unique_ptr<ParGridFunction> AvP;
    std::unique_ptr<ParGridFunction> AvB;


    
    double GetTotalPsi() const { return gtPsi; }
    double GetTotalPse() const { return gtPse; }
    int GetNE() const { return nE; }
    int GetNC() const { return nC; }
    int GetNV() const { return nV; }

    double L_w;

    Array<int> nbc_w_bdr;
    Array<int> ess_tdof_list_w;
    Array<int> ess_tdof_list_e;
    std::unique_ptr<ParMesh> pmesh;
    std::unique_ptr<ParMesh> pmesh0;







private:
    // Member functions
    void InitializeMesh();
    void CalculateElementVolume(int nE, const std::unique_ptr<ParMesh>& pmesh, Vector& EVol);
    void ReadGlobalDistanceFunction(const std::unique_ptr<FiniteElementSpace>& fespace);
    void InterpolateDomainParameters(int nV, const std::shared_ptr<ParFiniteElementSpace>& fespace);
    void CalculateTotals(const std::unique_ptr<ParGridFunction>& GridFunction, int nV, int nE, int nC, const Vector& EVol, double& localTotal, double& globalTotal);
    void CalculateTotalPsi(int nV, int nE, int nC, const Vector& EVol);
    void CalculateTotalPse(int nV, int nE, int nC, const Vector& EVol);
    void CalculateTargetCurrent(double tPsi);
    void PrintMeshInfo();

    mfem::ParGridFunction test_f; // Declare test_f as a member variable


    // Member variables
    const char* mesh_file;
    const char* dsF_file;
    int order;
    double dh;
    double zeta;
    double eps;
    double rho;
    double Cr;
    double Onm;
    mfem::Vector EVol;

    // Computed values
    double gtPsi;
    double gtPse;
    double gTrgI;

    // Mesh and FE space
    std::unique_ptr<Mesh> gmesh;
    std::unique_ptr<FiniteElementSpace> gFespace;
    // std::unique_ptr<ParFiniteElementSpace> fespace;
    std::shared_ptr<mfem::ParFiniteElementSpace> fespace;
    std::unique_ptr<ParGridFunction> dsF;
    std::unique_ptr<GridFunction> gDsF;
    std::unique_ptr<ParGridFunction> psi;
    std::unique_ptr<ParGridFunction> pse;
    // std::unique_ptr<ParGridFunction> AvP;
    // std::unique_ptr<ParGridFunction> AvB;

    // Other variables
    int nV; // Number of vertices
    int nE; // Number of elements
    int nC; // Number of corner vertices
    double localTotal; // Total computed value
    double val; // Intermediate value
    double globalTotal; // Global sum
    double tPsi; // Target Psi value
    double tPse; // Target Pse value
    double trgI; // Target current value
};

#endif // MESH_HANDLER_HPP
