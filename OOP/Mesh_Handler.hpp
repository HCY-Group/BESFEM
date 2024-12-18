#ifndef MESH_HANDLER_HPP
#define MESH_HANDLER_HPP

#include "Constants.hpp"
#include "mfem.hpp"
#include <memory>

using namespace mfem;

extern double gTrgI;


class MeshHandler {

public:
    MeshHandler();

    void LoadMesh();
    void Save();
    void SetupBoundaryConditions(mfem::ParMesh *pmesh, mfem::ParFiniteElementSpace *fespace);

    std::unique_ptr<mfem::ParGridFunction> psi; ///< Solid phase potential
    std::unique_ptr<mfem::ParGridFunction> pse; ///< Electrolyte phase potential
    std::unique_ptr<mfem::ParGridFunction> AvP; ///< Particle surface area
    std::unique_ptr<mfem::ParGridFunction> AvB; ///< Boundary surface area

    int nV; // Number of vertices
    int nE; // Number of elements
    int nC; // Number of corner vertices

    double L_w;

    mfem::Array<int> nbc_w_bdr;
    mfem::Array<int> ess_tdof_list_w;
    mfem::Array<int> ess_tdof_list_e;

    mfem::Array<int> dbc_w_bdr;
    mfem::Array<int> dbc_e_bdr;

    std::unique_ptr<mfem::ParMesh> pmesh;

    mfem::Vector EVol;

    double gtPsi;
    double gtPse;


private:
    // Member functions
    void InitializeMesh();
    void CalculateElementVolume(const std::unique_ptr<ParMesh>& pmesh, Vector& EVol);
    void ReadGlobalDistanceFunction(const std::unique_ptr<FiniteElementSpace>& fespace);
    void InitializeGridFunctions(const std::shared_ptr<mfem::ParFiniteElementSpace>& fespace);
    void InterpolateDomainParameters(const std::shared_ptr<ParFiniteElementSpace>& fespace);
    void CalculateTotals(const mfem::ParGridFunction& grid_function, const mfem::Vector& element_volumes, double& local_total, double& global_total);
    void CalculateTotalPhaseField(const mfem::ParGridFunction& grid_function, double& total, double& global_total);
    void CalculatePhasePotentialsAndTargetCurrent();
    void CalculateTargetCurrent(double total_psi);
    void PrintMeshInfo();

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
    double localTotal; // Total computed value
    double val; // Intermediate value
    double globalTotal; // Global sum
    double tPsi; // Target Psi value
    double tPse; // Target Pse value
    double trgI; // Target current value

    // Mesh and FE space
    std::unique_ptr<mfem::Mesh> gmesh;
    std::unique_ptr<mfem::FiniteElementSpace> gFespace;
    std::shared_ptr<mfem::ParFiniteElementSpace> fespace;
    std::unique_ptr<mfem::ParMesh> pmesh0;
    std::unique_ptr<mfem::ParGridFunction> dsF;
    std::unique_ptr<mfem::GridFunction> gDsF;

};

#endif // MESH_HANDLER_HPP
