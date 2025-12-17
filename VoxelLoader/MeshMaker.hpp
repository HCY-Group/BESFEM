
#include "mfem.hpp"

//using namespace std;
//using namespace mfem;

class MeshMaker {
public:
	MeshMaker(std::vector<std::vector<std::vector<int>>> tiffdata); //constructor
	void MakeGlobalMesh();
	void MakeParallelMesh();
	void AssignGlobalValues();
	void Make_H1_FESpace(int order=1);
	void Make_H1_FESpace_Parallel(int order = 1);
	void Make_DG_FESpace_Parallel(int order = 1);
	void TestGetFE();
	
	mfem::Mesh* GetGlobalMesh() {return gmesh;}
	mfem::ParMesh* GetParallelMesh() {return pmesh;}
	mfem::FiniteElementSpace* GetGlobalFESpace() {return gFespace;}
	mfem::ParFiniteElementSpace* GetParallelFESpace() {return fespace;}
	mfem::ParFiniteElementSpace* GetParallelFESpace_DG() {return fespace_dg;}
	mfem::ParFiniteElementSpace* GetParallelFESpace_DGdim() {return dimfespace_dg;}

private:
	std::vector<std::vector<std::vector<int>>> data;
	mfem::Mesh *gmesh = nullptr;
	mfem::ParMesh *pmesh = nullptr;
	mfem::H1_FECollection *gFec = nullptr;
	mfem::H1_FECollection *fec = nullptr;
	mfem::DG_FECollection *fec_dg = nullptr;
	//GridFunction *gVox;
	mfem::FiniteElementSpace *gFespace = nullptr;
	mfem::ParFiniteElementSpace *fespace = nullptr;
	mfem::ParFiniteElementSpace *fespace_dg = nullptr;
	mfem::ParFiniteElementSpace *dimfespace_dg = nullptr;

};
