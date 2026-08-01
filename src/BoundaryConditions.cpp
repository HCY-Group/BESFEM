#include "../include/BESFEM_All.hpp"
#include "mfem.hpp"
#include <tiffio.h>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <limits>
#include <cmath>
#include <mpi.h>
#include <iomanip>
#include <sstream>

using namespace std;
using sim::CellMode;
using sim::Electrode;

BoundaryConditions::BoundaryConditions(Initialize_Geometry &geo, Domain_Parameters &para)
    : geometry(geo), domain_parameters(para), parallelMesh(*geo.parallelMesh), globalMesh(*geo.globalMesh), parfespace(*geo.parfespace),
      E_L2G(geo.E_L2G) {}

BoundaryConditions::~BoundaryConditions() {}

void BoundaryConditions::SetupBoundaryConditions(CellMode mode, Electrode electrode)
{

    myid = mfem::Mpi::WorldRank();

    int dim = parallelMesh.Dimension();

    if (dim == 3)
    {

        if (mfem::Mpi::WorldRank() == 0) {std::cout << "Setting up boundary conditions for 3D mesh" << std::endl;}

        if (mode == CellMode::HALF && electrode == Electrode::ANODE)
        {
            if (mfem::Mpi::WorldRank() == 0) {std::cout << "Setting up boundary conditions for Half Cell: ANODE" << std::endl;}

            // East Neumann Boundary Condition
            nbc_e_bdr.SetSize(parallelMesh.bdr_attributes.Max());
            nbc_e_bdr = 0;
            nbc_e_bdr[2] = 1;

            // East Dirichlet Boundary Condition
            dbc_e_bdr.SetSize(parallelMesh.bdr_attributes.Max());
            dbc_e_bdr = 0;
            dbc_e_bdr[2] = 1;

            // West Dirichlet Boundary Condition
            dbc_w_bdr.SetSize(parallelMesh.bdr_attributes.Max());
            dbc_w_bdr = 0;
            dbc_w_bdr[4] = 1;

            ess_tdof_list_w.SetSize(0);
            parfespace.GetEssentialTrueDofs(dbc_w_bdr, ess_tdof_list_w);

            ess_tdof_list_e.SetSize(0);
            parfespace.GetEssentialTrueDofs(dbc_e_bdr, ess_tdof_list_e);

            nbc_bdr.SetSize(parallelMesh.bdr_attributes.Max());
            dbc_bdr.SetSize(parallelMesh.bdr_attributes.Max());

            ess_tdof_list = ess_tdof_list_w;

            nbc_bdr = nbc_e_bdr;
            dbc_bdr = dbc_e_bdr;
        }
        else if (mode == CellMode::HALF && electrode == Electrode::CATHODE)
        {
            if (mfem::Mpi::WorldRank() == 0) {std::cout << "Setting up boundary conditions for Half Cell: CATHODE" << std::endl;}

            // West Neumann Boundary Condition
            nbc_w_bdr.SetSize(parallelMesh.bdr_attributes.Max());
            nbc_w_bdr = 0;
            nbc_w_bdr[4] = 1;

            // West Dirichlet Boundary Condition
            dbc_w_bdr.SetSize(parallelMesh.bdr_attributes.Max());
            dbc_w_bdr = 0;
            dbc_w_bdr[4] = 1;
            // dbc_w_bdr[1] = 1;
            // dbc_w_bdr[2] = 1;
            // dbc_w_bdr[3] = 1;

            // East Dirichlet Boundary Condition
            dbc_e_bdr.SetSize(parallelMesh.bdr_attributes.Max());
            dbc_e_bdr = 0;
            dbc_e_bdr[2] = 1;

            ess_tdof_list_e.SetSize(0);
            parfespace.GetEssentialTrueDofs(dbc_e_bdr, ess_tdof_list_e);

            ess_tdof_list = ess_tdof_list_e;

            nbc_bdr = nbc_w_bdr;
            dbc_bdr = dbc_w_bdr;
        }
        else
        {
            if (mfem::Mpi::WorldRank() == 0) {std::cout << "Setting up boundary conditions for Full Cell" << std::endl;}

            // West Dirichlet Boundary Condition
            dbc_w_bdr.SetSize(parallelMesh.bdr_attributes.Max());
            dbc_w_bdr = 0;
            dbc_w_bdr[4] = 1;

            ess_tdof_list_w.SetSize(0);
            parfespace.GetEssentialTrueDofs(dbc_w_bdr, ess_tdof_list_w);

            // East Dirichlet Boundary Condition
            dbc_e_bdr.SetSize(parallelMesh.bdr_attributes.Max());
            dbc_e_bdr = 0;
            dbc_e_bdr[2] = 1;

            ess_tdof_list_e.SetSize(0);
            parfespace.GetEssentialTrueDofs(dbc_e_bdr, ess_tdof_list_e);

            nbc_bdr.SetSize(parallelMesh.bdr_attributes.Max());
            nbc_bdr = 0;

            dbc_bdr.SetSize(parallelMesh.bdr_attributes.Max());
            dbc_bdr = 0;

            SetupPinnedDOF(parfespace);
        }
    }
    else if (dim == 2)
    {

        if (mfem::Mpi::WorldRank() == 0) {std::cout << "Setting up boundary conditions for 2D mesh" << std::endl;}

        if (mode == CellMode::HALF && electrode == Electrode::ANODE)
        {
            if (mfem::Mpi::WorldRank() == 0) {std::cout << "Setting up boundary conditions for Half Cell: ANODE" << std::endl;}
            // East Neumann Boundary Condition
            nbc_e_bdr.SetSize(parallelMesh.bdr_attributes.Max());
            nbc_e_bdr = 0;
            nbc_e_bdr[1] = 1;
            nbc_e_bdr[2] = 1;
            nbc_e_bdr[3] = 1;
            nbc_e_bdr[0] = 1;

            // East Dirichlet Boundary Condition
            dbc_e_bdr.SetSize(parallelMesh.bdr_attributes.Max());
            dbc_e_bdr = 0;
            dbc_e_bdr[1] = 1;
            dbc_e_bdr[2] = 1;
            dbc_e_bdr[3] = 1;
            dbc_e_bdr[0] = 1;

            // West Dirichlet Boundary Condition
            dbc_w_bdr.SetSize(parallelMesh.bdr_attributes.Max());
            dbc_w_bdr = 0;
            dbc_w_bdr[3] = 1;

            ess_tdof_list_w.SetSize(0);
            parfespace.GetEssentialTrueDofs(dbc_w_bdr, ess_tdof_list_w);

            

            ess_tdof_list = ess_tdof_list_w;

            nbc_bdr.SetSize(parallelMesh.bdr_attributes.Max());
            dbc_bdr.SetSize(parallelMesh.bdr_attributes.Max());

            nbc_bdr = nbc_e_bdr;
            dbc_bdr = dbc_e_bdr;
        }
        else if (mode == CellMode::HALF && electrode == Electrode::CATHODE)
        {
            if (mfem::Mpi::WorldRank() == 0) {std::cout << "Setting up boundary conditions for Half Cell: CATHODE" << std::endl;}

            // Neumann Boundary Condition - used for electrolye concentration
            nbc_w_bdr.SetSize(parallelMesh.bdr_attributes.Max());
            nbc_w_bdr = 0;
            nbc_w_bdr[3] = 1; 

            // West Dirichlet Boundary Condition - used for electrolyte potential 
            dbc_w_bdr.SetSize(parallelMesh.bdr_attributes.Max());
            dbc_w_bdr = 0;
            dbc_w_bdr[3] = 1; 

            ess_tdof_list_w.SetSize(0);
            parfespace.GetEssentialTrueDofs(dbc_w_bdr, ess_tdof_list_w);

            // East Dirichlet Boundary Condition - used for particle potential
            dbc_e_bdr.SetSize(parallelMesh.bdr_attributes.Max());
            dbc_e_bdr = 0;
            dbc_e_bdr[1] = 1;

            ess_tdof_list_e.SetSize(0);
            parfespace.GetEssentialTrueDofs(dbc_e_bdr, ess_tdof_list_e);

            if (myid == 0)
            std::cout << "ess_tdof_list_e size = " 
                    << ess_tdof_list_e.Size() << std::endl;

            ess_tdof_list = ess_tdof_list_e;

            nbc_bdr = nbc_w_bdr;
            dbc_bdr = dbc_w_bdr;
        }
        else
        {
            if (mfem::Mpi::WorldRank() == 0) {std::cout << "Setting up boundary conditions for Full Cell" << std::endl;}

            // West Dirichlet Boundary Condition
            dbc_w_bdr.SetSize(parallelMesh.bdr_attributes.Max());
            dbc_w_bdr = 0;
            dbc_w_bdr[3] = 1;

            ess_tdof_list_w.SetSize(0);
            parfespace.GetEssentialTrueDofs(dbc_w_bdr, ess_tdof_list_w);

            // East Dirichlet Boundary Condition
            dbc_e_bdr.SetSize(parallelMesh.bdr_attributes.Max());
            dbc_e_bdr = 0;
            dbc_e_bdr[1] = 1;

            ess_tdof_list_e.SetSize(0);
            parfespace.GetEssentialTrueDofs(dbc_e_bdr, ess_tdof_list_e);

            nbc_bdr.SetSize(parallelMesh.bdr_attributes.Max());
            nbc_bdr = 0;

            dbc_bdr.SetSize(parallelMesh.bdr_attributes.Max());
            dbc_bdr = 0;

            SetupPinnedDOF(parfespace);
        }
    }
}

void BoundaryConditions::SetupPinnedDOF(mfem::ParFiniteElementSpace &fespace)
{
    myid = mfem::Mpi::WorldRank();

    double cand_dist = std::numeric_limits<double>::infinity();
    int local_lVpp = SelectCenterPin(0.99, cand_dist);

    struct
    {
        double dist;
        int rank;
    } local_candidate, global_candidate;

    local_candidate.dist = (local_lVpp >= 0) ? cand_dist : std::numeric_limits<double>::infinity();
    local_candidate.rank = myid;

    MPI_Allreduce(&local_candidate,
                  &global_candidate,
                  1,
                  MPI_DOUBLE_INT,
                  MPI_MINLOC,
                  MPI_COMM_WORLD);

    pin = false;
    rkpp = global_candidate.rank;

    ess_tdof_marker.SetSize(fespace.GetTrueVSize());
    ess_tdof_marker = 0;

    if (!std::isfinite(global_candidate.dist))
    {
        if (myid == 0)
            std::cerr << "[SetupPinnedDOF] Failed to elect a pinned vertex.\n";

        ess_tdof_list.SetSize(0);
        return;
    }

    if (myid == global_candidate.rank)
    {
        MFEM_VERIFY(local_lVpp >= 0,
                    "The elected rank does not have a valid local pinned vertex.");

        mfem::Array<int> vdofs;
        fespace.GetVertexVDofs(local_lVpp, vdofs);

        MFEM_VERIFY(vdofs.Size() > 0,
                    "The selected pinned vertex has no finite element DOFs.");

        int ldof = vdofs[0];
        if (ldof < 0)
            ldof = -1 - ldof;

        int ltdof = fespace.GetLocalTDofNumber(ldof);

        MFEM_VERIFY(ltdof >= 0,
                    "The selected pinned DOF is not owned by the elected rank.");

        ess_tdof_marker[ltdof] = 1;
        pin = true;

        const double *X = parallelMesh.GetVertex(local_lVpp);
        int dim = parallelMesh.Dimension();

        std::cout << "Rank " << myid
                  << " pinning local vertex " << local_lVpp << std::endl;

        std::cout << "coordinates of pinned vertex: (" << X[0];
        for (int d = 1; d < dim; d++)
            std::cout << ", " << X[d];
        std::cout << ")\n";
    }

    fespace.MarkerToList(ess_tdof_marker, ess_tdof_list);
}

// for picking the center
int BoundaryConditions::SelectCenterPin(double threshold, double &out_dist2)
{
    myid = mfem::Mpi::WorldRank();

    if (!domain_parameters.pse)
    {
        if (myid == 0)
            std::cerr << "[SelectPin] pse not initialized!\n";
        out_dist2 = std::numeric_limits<double>::infinity();
        return -1;
    }

    const mfem::ParGridFunction &pse_field = *domain_parameters.pse;

    int nE  = parallelMesh.GetNE();
    int dim = parallelMesh.Dimension();

    mfem::Vector minv(dim), maxv(dim);
    globalMesh.GetBoundingBox(minv, maxv);

    mfem::Vector center(dim);
    for (int d = 0; d < dim; d++)
        center[d] = 0.5 * (minv[d] + maxv[d]); // find the center

    mfem::Vector ec(dim);
    mfem::Array<int> VTX;
    mfem::Array<int> vdofs;

    double best_dist2 = std::numeric_limits<double>::infinity(); // first make the best distance positive inf
    int best_lVpp = -1; // first make the best local vertex -1

    const double z_tol = 1e-12; // tolerance for comparing z to min/max

    // loop over elements and find what is inside electrolyte
    for (int ei = 0; ei < nE; ei++)
    {
        parallelMesh.GetElementVertices(ei, VTX);
        int nV = VTX.Size();

        bool inside = false;

        for (int vi = 0; vi < nV; vi++)
        {
            parfespace.GetVertexDofs(VTX[vi], vdofs);

            if (vdofs.Size() == 0)
                continue;

            int ldof = vdofs[0];
            if (ldof < 0)
                ldof = -1 - ldof;

            if (ldof < 0 || ldof >= pse_field.Size())
                continue;

            double val = pse_field(ldof);
            if (val >= threshold)
            {
                inside = true;
                break;
            }
        }

        if (!inside) continue;

        // Use the refined parallel mesh directly. The serial global mesh is not
        // guaranteed to have the same element numbering after AMR.
        parallelMesh.GetElementCenter(ei, ec);

        if (dim == 3)
        {
            double zc = ec[2];
            double zmin = minv[2];
            double zmax = maxv[2];

            if (std::abs(zc - zmin) < z_tol || std::abs(zc - zmax) < z_tol)
                continue;
        }

        double dx = ec[0] - center[0]; // distance in x from global center
        double dy = (dim > 1) ? (ec[1] - center[1]) : 0.0; // distance in y from global center
        double dist2 = dx*dx + dy*dy;

        if (dist2 >= best_dist2)
            continue;

        // Find a vertex whose true DOF is owned by this rank. This prevents
        // electing a shared or constrained vertex that this rank cannot pin.
        int candidate_lVpp = -1;

        for (int vi = 0; vi < nV; vi++)
        {
            const double *X = parallelMesh.GetVertex(VTX[vi]);

            if (dim == 3)
            {
                double z = X[2];
                if (std::abs(z - minv[2]) < z_tol || std::abs(z - maxv[2]) < z_tol)
                    continue;
            }

            parfespace.GetVertexVDofs(VTX[vi], vdofs);

            if (vdofs.Size() == 0)
                continue;

            int ldof = vdofs[0];
            if (ldof < 0)
                ldof = -1 - ldof;

            int ltdof = parfespace.GetLocalTDofNumber(ldof);
            if (ltdof >= 0)
            {
                candidate_lVpp = VTX[vi];
                break;
            }
        }

        if (candidate_lVpp >= 0)
        {
            best_dist2 = dist2;
            best_lVpp = candidate_lVpp;
        }
    }

    if (best_lVpp < 0)
    {
        out_dist2 = std::numeric_limits<double>::infinity();
        return -1;
    }

    out_dist2 = best_dist2;
    return best_lVpp; // LOCAL vertex index candidate on this rank
}


int BoundaryConditions::SelectFirstPin(double threshold, double &out_dist2)
{
    myid = mfem::Mpi::WorldRank();

    if (!domain_parameters.pse)
    {
        if (myid == 0)
            std::cerr << "[SelectFirstPin] pse not initialized!\n";
        out_dist2 = std::numeric_limits<double>::infinity();
        return -1;
    }

    const mfem::ParGridFunction &pse_field = *domain_parameters.pse;

    int nE  = parallelMesh.GetNE();
    int dim = parallelMesh.Dimension();

    // ------------------------------
    // Bounding box for boundary check
    // ------------------------------
    mfem::Vector minv(dim), maxv(dim);
    globalMesh.GetBoundingBox(minv, maxv);

    const double tol = 1e-12;
    // const double tol = Constants::dh;

    mfem::Vector ec(dim);
    mfem::Array<int> VTX, vdofs;

    // ------------------------------
    // Loop over elements
    // ------------------------------
    for (int ei = 0; ei < nE; ei++)
    {
        parallelMesh.GetElementVertices(ei, VTX);
        int nV = VTX.Size();

        // --- Check pse to see if inside electrolyte ---
        bool inside = false;

        for (int vi = 0; vi < nV; vi++)
        {
            parfespace.GetVertexDofs(VTX[vi], vdofs);

            if (vdofs.Size() == 0)
                continue;

            int ldof = vdofs[0];
            if (ldof < 0)
                ldof = -1 - ldof;

            if (ldof < 0 || ldof >= pse_field.Size())
                continue;

            double val = pse_field(ldof);

            if (val >= threshold)
            {
                inside = true;
                break;
            }
        }

        if (!inside) continue;

        // --- Element center from the refined parallel mesh ---
        parallelMesh.GetElementCenter(ei, ec);

        double xc = ec[0];
        double yc = (dim > 1 ? ec[1] : 0.0);

        // Compute local element cell sizes dh_x, dh_y, dh_z from its vertices
        double ex_min = 1e300, ex_max = -1e300;
        double ey_min = 1e300, ey_max = -1e300;
        double ez_min = 1e300, ez_max = -1e300;

        for (int vi = 0; vi < nV; vi++)
        {
            const double *Xv = parallelMesh.GetVertex(VTX[vi]);
            ex_min = std::min(ex_min, Xv[0]);
            ex_max = std::max(ex_max, Xv[0]);
            if (dim > 1)
            {
                ey_min = std::min(ey_min, Xv[1]);
                ey_max = std::max(ey_max, Xv[1]);
            }
            if (dim > 2)
            {
                ez_min = std::min(ez_min, Xv[2]);
                ez_max = std::max(ez_max, Xv[2]);
            }
        }

        double dh_x = ex_max - ex_min;
        double dh_y = (dim > 1 ? ey_max - ey_min : 1.0);
        double dh_z = (dim > 2 ? ez_max - ez_min : 1.0);

        // Reject if center is within four local element widths of a boundary
        bool boundary_check = ((xc - minv[0]) < 4.0 * dh_x) || ((maxv[0] - xc) < 4.0 * dh_x) ||
            (dim > 1 && ((yc - minv[1]) < 4.0 * dh_y || (maxv[1] - yc) < 4.0 * dh_y));

        if (dim > 2)
        {
            double zc = ec[2];
            boundary_check = boundary_check ||
                ((zc - minv[2]) < 4.0 * dh_z) || ((maxv[2] - zc) < 4.0 * dh_z);
        }

        if (boundary_check)
            continue;

        // ------------------------------
        // FIRST valid element found!
        // Now select the first interior vertex owned by this rank
        // ------------------------------
        for (int vi = 0; vi < nV; vi++)
        {
            const double *X = parallelMesh.GetVertex(VTX[vi]);
            double x = X[0];
            double y = (dim > 1 ? X[1] : 0.0);
            double z = (dim > 2 ? X[2] : 0.0);

            bool vtx_on_boundary =
                (std::abs(x - minv[0]) < tol) || (std::abs(x - maxv[0]) < tol) ||
                (dim > 1 && (std::abs(y - minv[1]) < tol || std::abs(y - maxv[1]) < tol)) ||
                (dim > 2 && (std::abs(z - minv[2]) < tol || std::abs(z - maxv[2]) < tol));

            if (vtx_on_boundary)
                continue;

            parfespace.GetVertexVDofs(VTX[vi], vdofs);

            if (vdofs.Size() == 0)
                continue;

            int ldof = vdofs[0];
            if (ldof < 0)
                ldof = -1 - ldof;

            int ltdof = parfespace.GetLocalTDofNumber(ldof);
            if (ltdof >= 0)
            {
                out_dist2 = 0.0; // For MPI comparability
                return VTX[vi]; // LOCAL vertex index
            }
        }

        // If all vertices lie on the boundary or are not owned, skip element
    }

    // No interior element found
    out_dist2 = std::numeric_limits<double>::infinity();
    return -1;
}


void BoundaryConditions::SaveBoundaryConditionFields()
{
    myid = mfem::Mpi::WorldRank();

    mfem::ParGridFunction west_bc(&parfespace);
    mfem::ParGridFunction east_bc(&parfespace);
    mfem::ParGridFunction pinned_bc(&parfespace);

    west_bc = 0.0;
    east_bc = 0.0;
    pinned_bc = 0.0;

    mfem::Vector west_true(parfespace.GetTrueVSize());
    mfem::Vector east_true(parfespace.GetTrueVSize());
    mfem::Vector pinned_true(parfespace.GetTrueVSize());

    west_true = 0.0;
    east_true = 0.0;
    pinned_true = 0.0;

    for (int i = 0; i < ess_tdof_list_w.Size(); i++)
    {
        int tdof = ess_tdof_list_w[i];

        if (tdof >= 0 && tdof < west_true.Size())
            west_true[tdof] = 1.0;
    }

    for (int i = 0; i < ess_tdof_list_e.Size(); i++)
    {
        int tdof = ess_tdof_list_e[i];

        if (tdof >= 0 && tdof < east_true.Size())
            east_true[tdof] = 1.0;
    }

    for (int i = 0; i < ess_tdof_list.Size(); i++)
    {
        int tdof = ess_tdof_list[i];

        if (tdof >= 0 && tdof < pinned_true.Size())
            pinned_true[tdof] = 1.0;
    }

    west_bc.SetFromTrueDofs(west_true);
    east_bc.SetFromTrueDofs(east_true);
    pinned_bc.SetFromTrueDofs(pinned_true);

    std::ofstream mesh_ofs;
    std::ofstream west_ofs;
    std::ofstream east_ofs;
    std::ofstream pinned_ofs;

    if (myid == 0)
    {
        mesh_ofs.open("boundary_mesh");
        west_ofs.open("west_bc");
        east_ofs.open("east_bc");
        pinned_ofs.open("pinned_bc");

        mesh_ofs.precision(8);
        west_ofs.precision(8);
        east_ofs.precision(8);
        pinned_ofs.precision(8);
    }

    // Every rank must enter these collective calls.
    parallelMesh.PrintAsOne(mesh_ofs);
    west_bc.SaveAsOne(west_ofs);
    east_bc.SaveAsOne(east_ofs);
    pinned_bc.SaveAsOne(pinned_ofs);

    if (myid == 0)
    {
        mesh_ofs.close();
        west_ofs.close();
        east_ofs.close();
        pinned_ofs.close();

        std::cout << "Saved boundary_mesh, west_bc, east_bc, and pinned_bc."
                  << std::endl;
    }
}