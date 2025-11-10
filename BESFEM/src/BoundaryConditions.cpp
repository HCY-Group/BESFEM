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

using namespace std;
using sim::CellMode;
using sim::Electrode;



BoundaryConditions::BoundaryConditions(Initialize_Geometry& geo, Domain_Parameters &para)
    : geometry(geo), domain_parameters(para), parallelMesh(geo.parallelMesh.get()), globalMesh(geo.globalMesh.get()), parfespace(geo.parfespace),
    E_L2G(geo.E_L2G) {}

void BoundaryConditions::SetupBoundaryConditions(CellMode mode, Electrode electrode) {

    myid = mfem::Mpi::WorldRank();

    int dim = parallelMesh->Dimension();

    if (dim == 3) {

    std::cout << "Setting up boundary conditions for 3D mesh" << std::endl;

    if (mode == CellMode::HALF && electrode == Electrode::ANODE) {
        std::cout << "Setting up boundary conditions for Half Cell: ANODE" << std::endl;

        // East Neumann Boundary Condition
        nbc_e_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        nbc_e_bdr = 0;
        nbc_e_bdr[2] = 1; 

        // East Dirichlet Boundary Condition
        dbc_e_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        dbc_e_bdr = 0;
        dbc_e_bdr[2] = 1;

        // West Dirichlet Boundary Condition
        dbc_w_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        dbc_w_bdr = 0;
        dbc_w_bdr[0] = 1;

        ess_tdof_list_w.SetSize(0);
        parfespace->GetEssentialTrueDofs(dbc_w_bdr, ess_tdof_list_w);

        nbc_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        dbc_bdr.SetSize(parallelMesh->bdr_attributes.Max());

        nbc_bdr = nbc_e_bdr;
        dbc_bdr = dbc_e_bdr;


    } else if (mode == CellMode::HALF && electrode == Electrode::CATHODE) {
        std::cout << "Setting up boundary conditions for Half Cell: CATHODE" << std::endl;

        // West Neumann Boundary Condition
        nbc_w_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        nbc_w_bdr = 0;
        nbc_w_bdr[0] = 1; 

        // West Dirichlet Boundary Condition
        dbc_w_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        dbc_w_bdr = 0;
        dbc_w_bdr[0] = 1;

        // East Dirichlet Boundary Condition 
        dbc_e_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        dbc_e_bdr = 0;
        dbc_e_bdr[2] = 1;

        ess_tdof_list_e.SetSize(0);
        parfespace->GetEssentialTrueDofs(dbc_e_bdr, ess_tdof_list_e);

        nbc_bdr = nbc_w_bdr;
        dbc_bdr = dbc_w_bdr;

    } else {
        std::cout << "Setting up boundary conditions for Full Cell" << std::endl;

        // West Dirichlet Boundary Condition
        dbc_w_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        dbc_w_bdr = 0;
        dbc_w_bdr[0] = 1;

        ess_tdof_list_w.SetSize(0);
        parfespace->GetEssentialTrueDofs(dbc_w_bdr, ess_tdof_list_w);

        // East Dirichlet Boundary Condition
        dbc_e_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        dbc_e_bdr = 0;
        dbc_e_bdr[2] = 1;

        ess_tdof_list_e.SetSize(0);
        parfespace->GetEssentialTrueDofs(dbc_e_bdr, ess_tdof_list_e);

        nbc_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        nbc_bdr = 0;

        dbc_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        dbc_bdr = 0;

        SetupPinnedDOF(*parfespace);

    }


    } else if (dim == 2) {

    std::cout << "Setting up boundary conditions for 2D mesh" << std::endl;

    if (mode == CellMode::HALF && electrode == Electrode::ANODE) {
        std::cout << "Setting up boundary conditions for Half Cell: ANODE" << std::endl;

        // East Neumann Boundary Condition
        nbc_e_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        nbc_e_bdr = 0;
        nbc_e_bdr[2] = 1; 

        // East Dirichlet Boundary Condition
        dbc_e_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        dbc_e_bdr = 0;
        dbc_e_bdr[2] = 1;

        // West Dirichlet Boundary Condition
        dbc_w_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        dbc_w_bdr = 0;
        dbc_w_bdr[0] = 1;

        ess_tdof_list_w.SetSize(0);
        parfespace->GetEssentialTrueDofs(dbc_w_bdr, ess_tdof_list_w);

        nbc_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        dbc_bdr.SetSize(parallelMesh->bdr_attributes.Max());

        nbc_bdr = nbc_e_bdr;
        dbc_bdr = dbc_e_bdr;


    } else if (mode == CellMode::HALF && electrode == Electrode::CATHODE) {
        std::cout << "Setting up boundary conditions for Half Cell: CATHODE" << std::endl;

        // West Neumann Boundary Condition
        nbc_w_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        nbc_w_bdr = 0;
        nbc_w_bdr[0] = 1; 

        // West Dirichlet Boundary Condition
        dbc_w_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        dbc_w_bdr = 0;
        dbc_w_bdr[0] = 1;

        // East Dirichlet Boundary Condition 
        dbc_e_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        dbc_e_bdr = 0;
        dbc_e_bdr[2] = 1;

        ess_tdof_list_e.SetSize(0);
        parfespace->GetEssentialTrueDofs(dbc_e_bdr, ess_tdof_list_e);

        nbc_bdr = nbc_w_bdr;
        dbc_bdr = dbc_w_bdr;

    } else {
        std::cout << "Setting up boundary conditions for Full Cell" << std::endl;

        // West Dirichlet Boundary Condition
        dbc_w_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        dbc_w_bdr = 0;
        dbc_w_bdr[0] = 1;

        ess_tdof_list_w.SetSize(0);
        parfespace->GetEssentialTrueDofs(dbc_w_bdr, ess_tdof_list_w);

        // East Dirichlet Boundary Condition
        dbc_e_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        dbc_e_bdr = 0;
        dbc_e_bdr[2] = 1;

        ess_tdof_list_e.SetSize(0);
        parfespace->GetEssentialTrueDofs(dbc_e_bdr, ess_tdof_list_e);

        nbc_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        nbc_bdr = 0;

        dbc_bdr.SetSize(parallelMesh->bdr_attributes.Max());
        dbc_bdr = 0;

        SetupPinnedDOF(*parfespace);



    }

}

}

void BoundaryConditions::SetupPinnedDOF(mfem::ParFiniteElementSpace &fespace)
{
    // std::cout << "Rank: " << myid << " Setting up pinned DOF" << std::endl;

    if (E_L2G.Size() == 0)
        parallelMesh->GetGlobalElementIndices(E_L2G);


    // Pick gVpp automatically from electrolyte region
    int gVpp = ListElectrolyteElementVertices(0.9);
    if (gVpp < 0)
    {
        if (myid == 0)
            std::cerr << "[SetupPinnedDOF] No valid vertex found in electrolyte region!\n";
        return;
    }

    // int gVpp = 1311;
	int lVpp = -10;
	rkpp = -20;
	pin = false;

	// Map local distance function from global one
	for (int ei = 0; ei < geometry.nE; ei++){
		gei = E_L2G[ei];
		globalMesh->GetElementVertices(gei,gVTX);
		parallelMesh->GetElementVertices(ei,VTX);
		
		// identify pin point by comparing to global index
		for (int vi = 0; vi < geometry.nC; vi++){
			if (gVTX[vi] == gVpp){
				lVpp = VTX[vi];
				rkpp = myid;
				pin = true; 	
			}				
		}	
	}

	// imposing pinnig condition to the rank that has the pin point
	if (pin){
		mfem::Array<int> vdofs;
		fespace.GetVertexVDofs(lVpp, vdofs);	

		// mfem::Array<int> ess_tdof_marker(fespace.GetTrueVSize());
        ess_tdof_marker.SetSize(fespace.GetTrueVSize());
		ess_tdof_marker = 0;		
		
		int ldof = vdofs[0];
		if (ldof < 0) ldof = -1 - ldof;
		int ltdof = fespace.GetLocalTDofNumber(ldof); // -1 if this proc doesn't own the t-dof
		if (ltdof >= 0) ess_tdof_marker[ltdof] = 1;

		fespace.MarkerToList(ess_tdof_marker, ess_tdof_listPinned);
	}

    // std::cout << "Rank: " << rkpp << " Pinning global vertex " << gVpp << " as local vertex " << lVpp << std::endl;
	

    if (myid == rkpp){
        std::cout << "Rank " << myid << " found global vertex " << gVpp << " as local vertex " << lVpp << std::endl;
    }


}

int BoundaryConditions::ListElectrolyteElementVertices(double threshold)
{
    myid = mfem::Mpi::WorldRank();

    if (!domain_parameters.pse)
    {
        if (myid == 0)
            std::cerr << "[ListElectrolyteElementVertices] pse not initialized!\n";
        return -1;
    }

    const mfem::ParGridFunction &pse_field = *domain_parameters.pse;

    int nE = parallelMesh->GetNE();
    int dim = parallelMesh->Dimension();
    int nV_per_elem = std::pow(2, dim);
    mfem::Array<int> VTX(nV_per_elem);

    for (int ei = 0; ei < nE; ei++)
    {
        parallelMesh->GetElementVertices(ei, VTX);

        bool inside = true;
        for (int vi = 0; vi < VTX.Size(); vi++)
        {
            if (pse_field(VTX[vi]) < threshold)
            {
                inside = false;
                break;
            }
        }

        if (inside)
        {
            // pick the first vertex of the first valid element
            int chosen_vertex = VTX[0];
            if (myid == 0)
                std::cout << "[Rank " << myid << "] Selected vertex " << chosen_vertex
                          << " from electrolyte element " << ei << std::endl;
            return chosen_vertex;  // global vertex index
        }
    }

    return -1; // nothing found
}