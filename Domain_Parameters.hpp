#ifndef DOMAIN_PARAMETERS_HPP
#define DOMAIN_PARAMETERS_HPP

#include "mfem.hpp"
#include <memory>
#include <string>
#include <vector>

using namespace std;

extern double gTrgI;

class Domain_Parameters {

public:

    Domain_Parameters(Initialize_Geometry &geo);

    Initialize_Geometry &geometry;
    virtual ~Domain_Parameters();

    void SetupDomainParameters(const std::shared_ptr<mfem::ParFiniteElementSpace>& fespace);


    std::unique_ptr<mfem::ParGridFunction> psi; ///< Solid phase potential
    std::unique_ptr<mfem::ParGridFunction> pse; ///< Electrolyte phase potential
    std::unique_ptr<mfem::ParGridFunction> AvP; ///< Particle surface area
    std::unique_ptr<mfem::ParGridFunction> AvB; ///< Boundary surface area

    double gtPsi; ///< Global total for Psi
    double gtPse; ///< Global total for Pse


private:

    void InitializeGridFunctions(const std::shared_ptr<mfem::ParFiniteElementSpace>& fespace);
    void InterpolateDomainParameters(const std::shared_ptr<mfem::ParFiniteElementSpace>& fespace);
    void CalculateTotals(const mfem::ParGridFunction& grid_function, const mfem::Vector& element_volumes, double& local_total, double& global_total);
    void CalculateTotalPhaseField(const mfem::ParGridFunction& grid_function, double& total, double& global_total);
    void CalculatePhasePotentialsAndTargetCurrent();
    void CalculateTargetCurrent(double total_psi);
    void PrintInfo();

    mfem::Vector EVol; ///< Element volumes

    int nV;
    int nE;
    int nC;
    
    std::unique_ptr<mfem::ParGridFunction> dsF; ///< distance function grid
    std::unique_ptr<mfem::ParMesh> pmesh;        // Parallel mesh


    double tPsi; ///< Target Psi value
    double tPse; ///< Target Pse value
    double trgI; ///< Target current


};




#endif //DOMAIN_PARAMETERS_HPP