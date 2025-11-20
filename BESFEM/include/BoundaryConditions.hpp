
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
        virtual ~BoundaryConditions();

        void SetupBoundaryConditions(sim::CellMode mode, sim::Electrode electrode);
        void SetupPinnedDOF(mfem::ParFiniteElementSpace &fespace);
        int SelectCenterPin(double threshold, double &out_dist2);
        int SelectFirstPin(double threshold, double &out_dist2);


        // int SelectPin(double threshold, double &out_z, double &out_dist2);
        int gei; ///< Global element index.
        mfem::Array<int> nbc_w_bdr; ///< West Neumann Boundary Conditions
        mfem::Array<int> nbc_e_bdr; ///< East Neumann Boundary Conditions
        
        mfem::Array<int> nbc_bdr; ///< Electrolyte Neumann Boundary Conditions
        mfem::Array<int> dbc_bdr; ///< Electrolyte Dirichlet Boundary Conditions
        mfem::Array<int> dbc_w_bdr; ///< West Dirichlet Boundary Conditions
        mfem::Array<int> dbc_e_bdr; ///< East Dirichlet Boundary Conditions
        mfem::Array<int> ess_tdof_list_w; ///< Total DOF West
        mfem::Array<int> ess_tdof_list_e; ///< Total DOF East
        mfem::Array<int> gVTX;    ///< Global vertex indices of corners.
        mfem::Array<int> VTX;     ///< Local vertex indices of corners.
        // std::unique_ptr<mfem::Mesh> globalMesh;              ///< Global serial mesh.
        // std::shared_ptr<mfem::ParMesh> parallelMesh;
        // std::shared_ptr<mfem::ParFiniteElementSpace> parfespace; ///< Parallel finite element space.

        mfem::Mesh globalMesh;
        mfem::ParMesh parallelMesh;
        mfem::ParFiniteElementSpace parfespace;
        
        mfem::Array<HYPRE_BigInt> E_L2G; ///< Local to global element mapping.
        bool pin;
        mfem::Array<int> ess_tdof_list;
        mfem::Array<int> ess_tdof_marker; ///< East Dirichlet Boundary Conditions

        int myid; ///< MPI rank ID.
        int rkpp; ///< Rank owning the pinned DOF.
};
#endif // BOUNDARYCONDITIONS_HPP
