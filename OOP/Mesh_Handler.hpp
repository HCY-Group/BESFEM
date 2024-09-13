#ifndef MESH_HANDLER_HPP
#define MESH_HANDLER_HPP

#include "Constants.hpp"
#include <mfem.hpp>
#include <memory>

using namespace mfem;

class MeshHandler {
public:
    // Constructor
    MeshHandler();

    // Public member functions
    void LoadMesh();
    void Save();
    
private:
    // Member functions
    void InitializeMesh();
    void CalculateElementVolume(int nE, const std::unique_ptr<ParMesh>& pmesh, Vector& EVol);
    void ReadGlobalDistanceFunction(const std::unique_ptr<FiniteElementSpace>& fespace);
    void InterpolateDomainParameters(int nV, const std::unique_ptr<ParFiniteElementSpace>& fespace);
    void CalculateTotals(const std::unique_ptr<ParGridFunction>& GridFunction, int nV, int nE, int nC, const Vector& EVol, double& localTotal, double& globalTotal);
    void CalculateTotalPsi(int nV, int nE, int nC, const Vector& EVol);
    void CalculateTotalPse(int nV, int nE, int nC, const Vector& EVol);
    void CalculateTargetCurrent(double tPsi);
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

    // Computed values
    double gtPsi;
    double gtPse;
    double gTrgI;

    // Mesh and FE space
    std::unique_ptr<Mesh> gmesh;
    std::unique_ptr<ParMesh> pmesh;
    std::unique_ptr<FiniteElementSpace> gFespace;
    std::unique_ptr<ParFiniteElementSpace> fespace;
    std::unique_ptr<ParGridFunction> dsF;
    std::unique_ptr<GridFunction> gDsF;
    std::unique_ptr<ParGridFunction> psi;
    std::unique_ptr<ParGridFunction> pse;
    std::unique_ptr<ParGridFunction> AvP;
    std::unique_ptr<ParGridFunction> AvB;

    // Other variables
    int nV; // Number of vertices
    int nE; // Number of elements
    int nC; // Number of corner vertices
    double L_w; // West boundary size
    double localTotal; // Total computed value
    double val; // Intermediate value
    double globalTotal; // Global sum
    double tPsi; // Target Psi value
    double tPse; // Target Pse value
    double trgI; // Target current value
};

#endif // MESH_HANDLER_HPP
