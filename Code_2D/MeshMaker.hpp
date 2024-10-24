
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
	void TestGetFE();
	
	mfem::Mesh* GetGlobalMesh() {return gmesh;}
	mfem::ParMesh* GetParallelMesh() {return pmesh;}
	mfem::FiniteElementSpace* GetGlobalFESpace() {return gFespace;}
	mfem::ParFiniteElementSpace* GetParallelFESpace() {return fespace;}

private:
	std::vector<std::vector<std::vector<int>>> data;
	mfem::Mesh *gmesh = nullptr;
	mfem::ParMesh *pmesh = nullptr;
	mfem::H1_FECollection *gFec = nullptr;
	mfem::H1_FECollection *fec = nullptr;
	//GridFunction *gVox;
	mfem::FiniteElementSpace *gFespace = nullptr;
	mfem::ParFiniteElementSpace *fespace = nullptr;
};
