// PotentialBase.hpp
#pragma once
#include "mfem.hpp"
#include "Initialize_Geometry.hpp"
#include "Domain_Parameters.hpp"
#include "BoundaryConditions.hpp"
#include "Utils.hpp"
#include "FEMOperators.hpp"
#include <memory>

class PotentialBase {
public:
    
    PotentialBase(Initialize_Geometry &geo, Domain_Parameters &para);

    virtual ~PotentialBase() = default;

    virtual void SetupField(mfem::ParGridFunction &ph, double initial_value, mfem::ParGridFunction &psx) = 0;
    virtual void AssembleSystem(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &potential) = 0;
    virtual void UpdatePotential(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror) = 0;


    int    nE = 0;                            ///< Number of elements.
    int    nC = 0;                            ///< Nodal count per element (corners).
    int    nV = 0;                            ///< Number of vertices (global).
    const mfem::Vector &EVol;                 ///< Element volumes (from domain parameters).
    Initialize_Geometry &geometry;            ///< Geometry/mesh context
    Domain_Parameters   &domain_parameters;   ///< Domain/physics parameters

    mfem::ParMesh *pmesh = nullptr;                                  ///< Parallel mesh

    std::shared_ptr<mfem::ParFiniteElementSpace> fespace;             ///< Parallel FE space
    mfem::ParGridFunction TmpF; ///< Temporary field for error calculations

    mfem::ParGridFunction px0; ///< Pointer to the potential grid function
    mfem::HypreParVector X0; ///< Solution vector



};













