
#include "mfem.hpp"

using namespace std;
using namespace mfem;

class MeshMaker {
public:
	MeshMaker(std::vector<std::vector<std::vector<int>>> tiffdata); //constructor
	void MakeGlobalMesh();
	void MakeParallelMesh();
	void AssignGlobalValues();
	void Make_H1_FESpace(int order=1);
	void Make_H1_FESpace_Parallel(int order = 1);
	
	Mesh* GetGlobalMesh() {return gmesh;}
	ParMesh* GetParallelMesh() {return pmesh;}
	FiniteElementSpace* GetGlobalFESpace() {return gFespace;}
	ParFiniteElementSpace* GetParallelFESpace() {return fespace;}

private:
	std::vector<std::vector<std::vector<int>>> data;
	Mesh *gmesh;
	ParMesh *pmesh;
	//GridFunction *gVox;
	FiniteElementSpace *gFespace;
	ParFiniteElementSpace *fespace;
};
