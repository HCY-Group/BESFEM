
#include "mfem.hpp"

using namespace std;
using namespace mfem;

class MeshMaker {
public:
	MeshMaker(std::vector<std::vector<std::vector<int>>> tiffdata); //constructor
	void MakeGlobalMesh();
	
	GridFunction* GetGlobalVox() {return &gVox;}
	Mesh* GetGlobalMesh() {return &gmesh;}
	//Mesh GetGlobalMesh() {return gmesh;}
	//Mesh* GetGlobalMesh() const {return gFespace.GetMesh();}

private:
	std::vector<std::vector<std::vector<int>>> data;
	Mesh gmesh;
	GridFunction gVox;
	FiniteElementSpace gFespace;
};
