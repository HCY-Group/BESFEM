#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include "mpi.h"
#include "TimeDepOpers.hpp"
#include "readtiff.h"

using namespace std;
using namespace mfem;





int main(int argc, char *argv[])
{
	// Initialize MPI and HYPRE
	Mpi::Init(argc, argv);
	Hypre::Init();
	
	// ======================================
	// READ IN TIFF FILE
	// ======================================
	cout << "LOADING IN TIFF" << endl;
	const char* tiffname ="II_1_bin.tif";
	Constraints args;
	args.Depth_begin = 0;	//only read in one slice for 2D data
	args.Depth_end = 1;	//only read in one slice for 2D data
	TIFFReader reader(tiffname,args);
	reader.readinfo();
	std::vector<std::vector<std::vector<int>>> tiffdata;
	tiffdata = reader.getImageData();
	
	cout << "tiff size: "            << tiffdata.size()         << endl;	
	cout << "tiff[0] size: "         << tiffdata[0].size()      << endl;
	cout << "tiff[0][0] size: "      << tiffdata[0][0].size()   << endl;
	cout << "tiffdata[0][0][0] = "   << tiffdata[0][0][0]       << endl;
	/*for (int i = 0; i < tiffdata.size(); i++) {
		for (int j = 0; j < tiffdata[i].size(); j++) {
			for (int k = 0; k < tiffdata[i][j].size(); k++) {
				cout << tiffdata[i][j][k] << endl;
			}
		}
	}
	*/


	// ======================================
	// SERIAL MESH AND GLOBAL DATA
	// ======================================
	cout << "MAKING MESH" << endl;
	int nz = tiffdata.size();
	int ny = tiffdata[0].size();
	int nx = tiffdata[0][0].size();
	double sx = nx;  //make dx = 1
	double sy = ny;  //make dy = 1
	bool generate_edges = false;
	bool sfc_ordering = false;
	Mesh gmesh = Mesh::MakeCartesian2D(nx-1, ny-1, Element::QUADRILATERAL, generate_edges, sx, sy, sfc_ordering);
	gmesh.EnsureNCMesh(true);

	// Create global FE space for Voxel Data
	int order = 1;
	H1_FECollection gFec(order, gmesh.Dimension());
	FiniteElementSpace gFespace(&gmesh, &gFec);

	// global grid function for voxel data
	cout << "Defining Voxel GridFunction" << endl;
	GridFunction gVox(&gFespace);
	for(int j=0; j<ny; j++){
		for(int i=0; i<nx; i++){
			int idx = i+ny*j;
			//cout << "idx = " << i+ny*j << endl;
			//cout << "gVox[idx] = " << gVox[idx] << endl;
			//cout << "tiffdata[i][j][0] = " << tiffdata[i][j][0] << endl;
			gVox[idx] = tiffdata[0][j][i];
		}
	}
	
	// double check our read-in: output gVox to Paraview
	cout << "PRINTING gVox" << endl;
	ParaViewDataCollection *pd = NULL;
	pd = new ParaViewDataCollection("gVoxelData", &gmesh);
	pd->RegisterField("gVox", &gVox);
	pd->SetLevelsOfDetail(order);
	pd->SetDataFormat(VTKFormat::BINARY);
	pd->SetHighOrderOutput(true);
	pd->SetCycle(0);
	pd->SetTime(0.0);
	pd->Save();
	delete pd;
		
	/*
	// glvis -m gmesh.mesh -g gvox.gf
	gVox.Save("gvox.gf");
	gmesh.Save("gmesh.mesh");
	*/

	/*
	char vishost[] = "localhost";
 	int  visport   = 19916;
      	socketstream sol_sock(vishost, visport);
      	sol_sock.precision(8);
      	sol_sock << "solution\n" << gmesh << gVox << flush;
	*/
	/*

	char vishost[] = "localhost";
 	int  visport   = 19916;
      	socketstream sol_sock(vishost, visport);
      	sol_sock << "parallel " << Mpi::WorldSize() << " " << Mpi::WorldRank() << "\n";
      	sol_sock.precision(8);
      	sol_sock << "solution\n" << gmesh << gVox << flush;
	*/











	// ======================================
	// LOCAL (PARALLEL) MESH AND DATA
	// ======================================
	ParMesh pmesh(MPI_COMM_WORLD, gmesh);
	int nV = pmesh.GetNV();			//number of vertices
	int nE = pmesh.GetNE();			//number of elements
	int nC = pow(2,pmesh.Dimension());	//number of corner vertices
	
	// Define finite element space
	H1_FECollection fec(order, pmesh.Dimension());
	ParFiniteElementSpace fespace(&pmesh, &fec);

	// Map local to global element indices
	Array<HYPRE_BigInt> E_L2G;
	pmesh.GetGlobalElementIndices(E_L2G);

	// Local GridFunction
	cout << "DEFINING LOCAL GRIDFUNCTION" << endl;
	ParGridFunction Vox(&fespace);
	
	Array<int> gVTX(nC);	//global indices of corner vertices
	Array<int> VTX(nC);	//local indices of corner vertices
	int gei;			//global element indices
	int ei;			//local element indices
	for (ei=0; ei<nE; ei++){
		gei = E_L2G[ei];

		gmesh.GetElementVertices(gei,gVTX);
		pmesh.GetElementVertices(ei,VTX);
	
		for (int vi = 0; vi<nC; vi++){
			Vox(VTX[vi]) = gVox(gVTX[vi]);
		}
	}
	


	// double check our read-in: output Vox to Paraview
	cout << "PRINTING OUT Vox" << endl;
	//ParaViewDataCollection *pd = NULL;
	pd = NULL;
	pd = new ParaViewDataCollection("pVoxelData", &pmesh);
	pd->RegisterField("Vox", &Vox);
	pd->SetLevelsOfDetail(order);
	pd->SetDataFormat(VTKFormat::BINARY);
	pd->SetHighOrderOutput(true);
	pd->SetCycle(0);
	pd->SetTime(0.0);
	pd->Save();
	delete pd;






	// ======================================
	// SMOOTH THE DATA USING ALLEN-CAHN
	// ======================================



	// ======================================
	// FIND DISTANCE FUNCTION BY REINITIALIZING LEVEL SET
	// ======================================


	// ======================================
	// FIND CONNECTIVITY
	// ======================================

	// Natural (Neumann) boundary conditions
	Array<int> boundary_dofs;

	return 0;

}
