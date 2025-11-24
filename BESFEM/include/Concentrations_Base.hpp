// ConcentrationBase.hpp
#pragma once
#include "mfem.hpp"
#include "Initialize_Geometry.hpp"
#include "Domain_Parameters.hpp"
#include "BoundaryConditions.hpp"
#include "Utils.hpp"
#include "FEMOperators.hpp"
#include <memory>

class ConcentrationBase {
public:
    
    ConcentrationBase(Initialize_Geometry &geo, Domain_Parameters &para);

    virtual ~ConcentrationBase() = default;

    virtual void SetupField(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx) = 0;

    virtual void UpdateConcentration(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx) = 0;

    double GetLithiation() const { return Xfr; } ///< Return current degree of lithiation.

    double Xfr = 0.0;
    double CeC = 0.0;
    double gCeC = 0.0;                                                ///< Global total salt
    mfem::Array<double> VtxVal;                                        ///< Vertex values buffer.


    int    nE = 0;                            ///< Number of elements.
    int    nC = 0;                            ///< Nodal count per element (corners).
    int    nV = 0;                            ///< Number of vertices (global).
    mfem::Vector EAvg;                        ///< Per-element averages (work buffer).
    const mfem::Vector &EVol;                 ///< Element volumes (from domain parameters).
    double CeAvg = 0.0;                                               ///< Global average salt.
    double Ce0   = 0.001;                                             ///< Initial salt reference.

    double gtPsi = 0.0;                                                ///< Global ψ_P normalization.
    double gtPse = 0.0;                                                ///< Global ψ_E normalization.
    mfem::Array<int> boundary_dofs; ///< Boundary degrees of freedom
    mfem::HypreParVector X1v; ///< Scratch vector

    Initialize_Geometry &geometry;            ///< Geometry/mesh context
    Domain_Parameters   &domain_parameters;   ///< Domain/physics parameters

    mfem::ParMesh *pmesh = nullptr;                                  ///< Parallel mesh
    mfem::Mesh    *gmesh = nullptr;                                   ///< Serial/global mesh

    std::shared_ptr<mfem::ParFiniteElementSpace> fespace;             ///< Parallel FE space
    std::unique_ptr<mfem::ParGridFunction> TmpF;                       ///< Temporary field for products.



};


