#ifndef BOUNDARYCONDITIONS_HPP
#define BOUNDARYCONDITIONS_HPP

#include "mfem.hpp"
#include <memory>
#include <string>
#include <vector>
#include "SimTypes.hpp"

class Initialize_Geometry;
class Domain_Parameters;

/**
 * @class BoundaryConditions
 * @brief Class for setting up boundary conditions in battery simulations.
*/
class BoundaryConditions {
    public: 

        BoundaryConditions(Initialize_Geometry &geo, Domain_Parameters &para);

        Initialize_Geometry &geometry;
        Domain_Parameters &domain_parameters;

        void SetupBoundaryConditions(sim::CellMode mode, sim::Electrode electrode);
        void SetupPinnedDOF(mfem::ParFiniteElementSpace &fespace);
        int ListElectrolyteElementVertices(double threshold);

        int gei; ///< Global element index.
        int ei;  ///< Local element index.

        mfem::Array<int> nbc_w_bdr; ///< West Neumann Boundary Conditions
        mfem::Array<int> nbc_s_bdr; ///< South Neumann Boundary Conditions
        mfem::Array<int> nbc_e_bdr; ///< East Neumann Boundary Conditions
        mfem::Array<int> nbc_n_bdr; ///< North Neumann Boundary Conditions
        
        mfem::Array<int> nbc_bdr; ///< Electrolyte Neumann Boundary Conditions
        mfem::Array<int> dbc_bdr; ///< Electrolyte Dirichlet Boundary Conditions

        mfem::Array<int> dbc_w_bdr; ///< West Dirichlet Boundary Conditions
        mfem::Array<int> dbc_e_bdr; ///< East Dirichlet Boundary Conditions

        mfem::Array<int> ess_tdof_list_w; ///< Total DOF West
        mfem::Array<int> ess_tdof_list_e; ///< Total DOF East

        mfem::Array<int> gVTX;    ///< Global vertex indices of corners.
        mfem::Array<int> VTX;     ///< Local vertex indices of corners.

        std::unique_ptr<mfem::Mesh> globalMesh;              ///< Global serial mesh.
        std::shared_ptr<mfem::ParMesh> parallelMesh;
        std::shared_ptr<mfem::FiniteElementSpace> feSpace;   ///< Serial finite element space.
        std::shared_ptr<mfem::FiniteElementSpace> globalfespace; ///< Global finite element space.
        std::shared_ptr<mfem::ParFiniteElementSpace> parfespace; ///< Parallel finite element space.
        std::shared_ptr<mfem::ParFiniteElementSpace> parfespace_dg; ///< Parallel DG finite element space.
        std::shared_ptr<mfem::ParFiniteElementSpace> pardimfespace_dg; ///< Parallel DG FE space (dim).
        mfem::Array<HYPRE_BigInt> E_L2G; ///< Local to global element mapping.

        double Onm; ///< Number of grid function entries.
        std::unique_ptr<mfem::GridFunction> gDsF; ///< Global distance function.
        std::unique_ptr<mfem::ParGridFunction> dsF; ///< Parallel distance function.

        std::unique_ptr<mfem::GridFunction>       gDsF_A;  ///< Global distance function for anode.
        std::unique_ptr<mfem::GridFunction>       gDsF_C;  ///< Global distance function for cathode.
        std::unique_ptr<mfem::ParGridFunction>    dsF_A;    ///< Parallel distance function for anode.
        std::unique_ptr<mfem::ParGridFunction>    dsF_C;    ///< Parallel distance function for cathode.
        
        std::unique_ptr<mfem::GridFunction> gVox; ///< Global voxel function.
        std::unique_ptr<mfem::ParGridFunction> Vox; ///< Parallel voxel function.

        std::vector<std::vector<std::vector<int>>> tiffData; ///< TIFF voxel data.
        std::unique_ptr<mfem::H1_FECollection> gfec; ///< Global H1 finite element collection.
        std::unique_ptr<mfem::H1_FECollection> pfec; ///< Parallel H1 finite element collection.
        std::unique_ptr<mfem::DG_FECollection> pfec_dg; ///< Parallel DG finite element collection.

        mfem::Array<int> ess_tdof_potE;        // size 0 on most ranks, size 1 on owner
        bool anchor_set = false;
        HYPRE_BigInt global_anchor_potE;  // the global true dof index for the anchor
        int anchor_owner_potE; 
        bool pin;

        mfem::Array<int> ess_tdof_listPinned;
        mfem::Array<int> boundary_dofs;
        mfem::Array<int> ess_tdof_marker; ///< East Dirichlet Boundary Conditions


        int myid; ///< MPI rank ID.
        int rkpp; ///< Rank owning the pinned DOF.

};

#endif // BOUNDARYCONDITIONS_HPP