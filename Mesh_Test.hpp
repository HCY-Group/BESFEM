#ifndef MESH_TEST_HPP
#define MESH_TEST_HPP

#include "MeshBase.hpp"
#include "mfem.hpp"
#include <memory>

class Mesh_Test : public MeshBase {
private:
    // std::unique_ptr<mfem::ParGridFunction> psi;    // Phase field grid function
    // std::unique_ptr<mfem::ParGridFunction> pse;    // Complementary phase field
    // mfem::Vector elementAverages;                 // Element average values
    // double targetCurrent = 0.0;                   // Target current
    // double totalPsi = 0.0, totalPse = 0.0;        // Total Psi and Pse values

public:
    Mesh_Test();
    ~Mesh_Test();

    void Initialize();
    // void SetupBoundaryConditions() override;
    // void CalculatePhaseFields();
    // void PrintResults() const;
};

#endif // MESH_TEST_HPP
