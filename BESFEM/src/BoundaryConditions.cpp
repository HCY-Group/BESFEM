#include "../include/BoundaryConditions.hpp"
#include "../include/Initialize_Geometry.hpp"
#include "../include/Domain_Parameters.hpp"
#include "../inputs/Constants.hpp"
#include "../include/readtiff.h"
#include "../include/SimTypes.hpp"
#include "mfem.hpp"
#include <tiffio.h>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <limits>
#include <cmath>
#include <mpi.h>

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

        std::cout << "Setting up boundary conditions for 3D mesh" << std::endl;

        if (mode == CellMode::HALF && electrode == Electrode::ANODE)
        {
            std::cout << "Setting up boundary conditions for Half Cell: ANODE" << std::endl;

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
            dbc_w_bdr[0] = 1;

            ess_tdof_list_w.SetSize(0);
            parfespace.GetEssentialTrueDofs(dbc_w_bdr, ess_tdof_list_w);

            nbc_bdr.SetSize(parallelMesh.bdr_attributes.Max());
            dbc_bdr.SetSize(parallelMesh.bdr_attributes.Max());

            ess_tdof_list = ess_tdof_list_w;

            nbc_bdr = nbc_e_bdr;
            dbc_bdr = dbc_e_bdr;
        }
        else if (mode == CellMode::HALF && electrode == Electrode::CATHODE)
        {
            std::cout << "Setting up boundary conditions for Half Cell: CATHODE" << std::endl;

            // West Neumann Boundary Condition
            nbc_w_bdr.SetSize(parallelMesh.bdr_attributes.Max());
            nbc_w_bdr = 0;
            nbc_w_bdr[0] = 1;

            // West Dirichlet Boundary Condition
            dbc_w_bdr.SetSize(parallelMesh.bdr_attributes.Max());
            dbc_w_bdr = 0;
            dbc_w_bdr[0] = 1;

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
            std::cout << "Setting up boundary conditions for Full Cell" << std::endl;

            // West Dirichlet Boundary Condition
            dbc_w_bdr.SetSize(parallelMesh.bdr_attributes.Max());
            dbc_w_bdr = 0;
            dbc_w_bdr[0] = 1;

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

        std::cout << "Setting up boundary conditions for 2D mesh" << std::endl;

        if (mode == CellMode::HALF && electrode == Electrode::ANODE)
        {
            std::cout << "Setting up boundary conditions for Half Cell: ANODE" << std::endl;

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
            dbc_w_bdr[0] = 1;

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
            std::cout << "Setting up boundary conditions for Half Cell: CATHODE" << std::endl;

            // West Neumann Boundary Condition
            nbc_w_bdr.SetSize(parallelMesh.bdr_attributes.Max());
            nbc_w_bdr = 0;
            nbc_w_bdr[0] = 1;

            // West Dirichlet Boundary Condition
            dbc_w_bdr.SetSize(parallelMesh.bdr_attributes.Max());
            dbc_w_bdr = 0;
            dbc_w_bdr[0] = 1;

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
            std::cout << "Setting up boundary conditions for Full Cell" << std::endl;

            // West Dirichlet Boundary Condition
            dbc_w_bdr.SetSize(parallelMesh.bdr_attributes.Max());
            dbc_w_bdr = 0;
            dbc_w_bdr[0] = 1;

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
}

void BoundaryConditions::SetupPinnedDOF(mfem::ParFiniteElementSpace &fespace)
{
    myid = mfem::Mpi::WorldRank();

    if (E_L2G.Size() == 0)
        parallelMesh.GetGlobalElementIndices(E_L2G);

    double cand_dist = 0.0;
    int local_gVpp = SelectFirstPin(0.99, cand_dist);

    // std::cout << "[SetupPinnedDOF] rank " << myid
    //           << " gVpp = "  << local_gVpp
    //           << ", element distance to center = " << cand_dist
    //           << std::endl;

    // const double *X = globalMesh.GetVertex(local_gVpp);
    // std::cout << "Rank: " << myid << " pinned coords: (" << X[0] << ", " << X[1] << ", " << X[2] << ")\n";

    const double posInf = std::numeric_limits<double>::infinity();

    if (local_gVpp < 0)
    {
        cand_dist = posInf;
    }

    double local_d = cand_dist;
    double best_d = 0.0;
    MPI_Allreduce(&local_d, &best_d, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD); // find the lowest global distance from center value

    const double d_tol = 1e-14;

    int filtered_vid = (std::abs(cand_dist - best_d) <= d_tol && local_gVpp >= 0) ? local_gVpp: std::numeric_limits<int>::max();

    int gVpp = 0;
    MPI_Allreduce(&filtered_vid, &gVpp, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD); // find the smallest gVpp after finding the smallest distance

    if (gVpp == std::numeric_limits<int>::max())
    {
        if (myid == 0)
            std::cerr << "[SetupPinnedDOF] Failed to elect a global pinned vertex.\n";
        return;
    }

    int lVpp = -1;
    pin = false;
    rkpp = -1;

    for (int ei = 0; ei < geometry.nE; ei++)
    {
        int gei = E_L2G[ei];
        globalMesh.GetElementVertices(gei, gVTX);
        parallelMesh.GetElementVertices(ei, VTX);

        for (int vi = 0; vi < geometry.nC; vi++)
        {
            if (gVTX[vi] == gVpp)
            {
                lVpp = VTX[vi]; // local vertex index on this rank
                pin = true;
                rkpp = myid;
            }
        }
    }

    if (pin)
    {
        mfem::Array<int> vdofs;
        fespace.GetVertexVDofs(lVpp, vdofs);

        ess_tdof_marker.SetSize(fespace.GetTrueVSize());
        ess_tdof_marker = 0;

        int ldof = vdofs[0];
        if (ldof < 0)
            ldof = -1 - ldof;

        int ltdof = fespace.GetLocalTDofNumber(ldof);
        if (ltdof >= 0)
            ess_tdof_marker[ltdof] = 1;

        fespace.MarkerToList(ess_tdof_marker, ess_tdof_list);

        std::cout << "Rank " << myid
                  << " pinning global vertex " << gVpp
                  << " as local vertex " << lVpp << std::endl;

        const double *X = globalMesh.GetVertex(gVpp);
        std::cout << "coordinates of pinned vertex: (" << X[0] << ", " << X[1] << ", " << X[2] << ")\n";
    }

}

// for picking the center
int BoundaryConditions::SelectCenterPin(double threshold, double &out_dist2)
{
    myid = mfem::Mpi::WorldRank();

    if (!domain_parameters.pse)
    {
        if (myid == 0)
            std::cerr << "[SelectPin] pse not initialized!\n";
        out_dist2 =  std::numeric_limits<double>::infinity();
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

    if (E_L2G.Size() == 0)
        parallelMesh.GetGlobalElementIndices(E_L2G);

    mfem::Vector ec(dim);
    mfem::Array<int> VTX;
    mfem::Array<int> vdofs;

    double best_dist2 =  std::numeric_limits<double>::infinity(); // first make the best distance positive inf
    int    best_ei    = -1; // first make the best element -1

    const double z_tol = 1e-12; // tolerance for comparing z to min/max

    // loop over vertices and find what is inside electrolyte
    for (int ei = 0; ei < nE; ei++)
    {
        parallelMesh.GetElementVertices(ei, VTX);
        int nV = VTX.Size();

        bool inside = false;

        if (dim == 2)
        {
            inside = false;
            for (int vi = 0; vi < nV; vi++)
            {
                if (pse_field(VTX[vi]) >= threshold)
                {
                    inside = true;
                    break;
                }
            }
        }
        else if (dim == 3)
        {
            inside = false;
            for (int vi = 0; vi < nV; vi++)
            {
                parfespace.GetVertexDofs(VTX[vi], vdofs);
                double val = pse_field(vdofs[0]);
                if (val >= threshold)
                {
                    inside = true;
                    break;
                }
            }
        }

        if (!inside) continue;

        int gei = E_L2G[ei];
        globalMesh.GetElementCenter(gei, ec);

        if (dim == 3)
        {
            double zc = ec[2];
            double zmin = minv[2];
            double zmax = maxv[2];

            if (std::abs(zc - zmin) < z_tol || std::abs(zc - zmax) < z_tol)
            {
                continue;
            }
        }

        double dx = ec[0] - center[0]; // distance in x from global center
        double dy = (dim > 1) ? (ec[1] - center[1]) : 0.0; // distance in y from global center
        double dist2 = dx*dx + dy*dy;

        // pick the element closest to the center in the xy-plane;
        // ties are automatically broken by "first one encountered"
        if (dist2 < best_dist2)
        {
            best_dist2 = dist2;
            best_ei    = ei;
        }

    }

    // all elements not the best one are -1
    if (best_ei < 0)
    {
        out_dist2 =  std::numeric_limits<double>::infinity();
        return -1;
    }

    // need to find which vertex to pin
    parallelMesh.GetElementVertices(best_ei, VTX);
    int nV = VTX.Size();
    int gei = E_L2G[best_ei];

    mfem::Array<int> gVTX(nV);
    globalMesh.GetElementVertices(gei, gVTX);

    int chosen_gv = gVTX[0]; // default

    if (dim == 3)
    {
        double zmin = minv[2];
        double zmax = maxv[2];
        const double z_tol = 1e-12;  // same tolerance as before

        chosen_gv = -1;

        // Prefer a vertex whose z is strictly between zmin and zmax
        for (int vi = 0; vi < nV; vi++)
        {
            const double *X = globalMesh.GetVertex(gVTX[vi]); // [x,y,z]
            double z = X[2];

            if (std::abs(z - zmin) >= z_tol && std::abs(z - zmax) >= z_tol)
            {
                chosen_gv = gVTX[vi];
                break; // first such vertex wins
            }
        }

        // if all vertices are on zmin or zmax (should be rare), just pick the first one so we still have something to pin
        if (chosen_gv < 0)
            chosen_gv = gVTX[0];
    }
    else
        {
            // 2D: z doesn't matter
            chosen_gv = gVTX[0];
        }

        out_dist2 = best_dist2;
        return chosen_gv;  // GLOBAL vertex index candidate on this rank

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

    if (E_L2G.Size() == 0)
        parallelMesh.GetGlobalElementIndices(E_L2G);

    mfem::Vector ec(dim);
    mfem::Array<int> VTX, gVTX, vdofs;

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
            double val = pse_field(vdofs[0]);

            if (val >= threshold)
            {
                inside = true;
                break;
            }
        }

        if (!inside) continue;

        // --- Element center in global coords ---
        int gei = E_L2G[ei];
        globalMesh.GetElementCenter(gei, ec);

        double xc = ec[0];
        double yc = (dim > 1 ? ec[1] : 0.0);
        double zc = (dim > 2 ? ec[2] : 0.0);

        // // --- Reject element touching ANY boundary ---
        // bool elem_on_boundary =

        //     // X boundaries
        //     (std::abs(xc - minv[0]) < tol) || (std::abs(xc - maxv[0]) < tol) ||

        //     // Y boundaries
        //     (dim > 1 && (std::abs(yc - minv[1]) < tol || std::abs(yc - maxv[1]) < tol)) ||

        //     // Z boundaries
        //     (dim > 2 && (std::abs(zc - minv[2]) < tol || std::abs(zc - maxv[2]) < tol));

        // if (elem_on_boundary)
        //     continue;

        // Compute local element cell sizes dh_x, dh_y, dh_z from its vertices
        double ex_min = 1e300, ex_max = -1e300;
        double ey_min = 1e300, ey_max = -1e300;
        double ez_min = 1e300, ez_max = -1e300;

        globalMesh.GetElementVertices(gei, gVTX);
        for (int vi = 0; vi < gVTX.Size(); vi++)
        {
            const double *Xv = globalMesh.GetVertex(gVTX[vi]);
            ex_min = std::min(ex_min, Xv[0]);
            ex_max = std::max(ex_max, Xv[0]);
            if (dim > 1) {
                ey_min = std::min(ey_min, Xv[1]);
                ey_max = std::max(ey_max, Xv[1]);
            }
            if (dim > 2) {
                ez_min = std::min(ez_min, Xv[2]);
                ez_max = std::max(ez_max, Xv[2]);
            }
        }

        double dh_x = ex_max - ex_min;
        double dh_y = (dim > 1 ? ey_max - ey_min : 1.0);
        double dh_z = (dim > 2 ? ez_max - ez_min : 1.0);

        // Reject if center within certain number of elements in x and y direction
        // if the element is too close to and x or y direction, then I keep looking for another element
        bool boundary_check = ((xc - minv[0]) < 4.0 * dh_x) || ((maxv[0] - xc) < 4.0 * dh_x) ||
            (dim > 1 && ((yc - minv[1]) < 4.0 * dh_y || (maxv[1] - yc) < 4.0 * dh_y));

        if (boundary_check) // if the element ends up not being too close, then do not continue and keep that for the pin
            continue;


        // ------------------------------
        // FIRST valid element found!
        // Now select the first interior vertex of this element
        // ------------------------------
        globalMesh.GetElementVertices(gei, gVTX);

        // pick vertex that is interior in ALL dims
        for (int vi = 0; vi < nV; vi++)
        {
            const double *X = globalMesh.GetVertex(gVTX[vi]);
            double x = X[0];
            double y = (dim > 1 ? X[1] : 0.0);
            double z = (dim > 2 ? X[2] : 0.0);

            bool vtx_on_boundary =
                (std::abs(x - minv[0]) < tol) || (std::abs(x - maxv[0]) < tol) ||
                (dim > 1 && (std::abs(y - minv[1]) < tol || std::abs(y - maxv[1]) < tol)) ||
                (dim > 2 && (std::abs(z - minv[2]) < tol || std::abs(z - maxv[2]) < tol));

            if (!vtx_on_boundary)
            {
                out_dist2 = 0.0;   // For MPI comparability
                return gVTX[vi];   // GLOBAL vertex index
            }
        }

        // If all vertices lie on boundary, skip element
        // (rare but possible with structured meshes)
    }

    // No interior element found
    out_dist2 = std::numeric_limits<double>::infinity();
    return -1;
}
