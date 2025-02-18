
CONFIG_MK = $(EBROOTMFEM)/share/mfem/config.mk

include $(CONFIG_MK)


LDFLAGS = -ltiff
SOURCECODE = ../code/Constants.cpp ../code/Initialize_Geometry.cpp TimeDepOpers.cpp readtiff.cpp MeshMaker.cpp VoxelSolver.cpp VoxelSolver_DG.cpp

MFEM_CXXFLAGS = -std=c++14 -O3 -DNDEBUG -O2 -ftree-vectorize -march=native -fno-math-errno -fopenmp

# Remove built-in rule
%: %.cpp

# Replace the default implicit rule for *.cpp files
%: %.cpp $(MFEM_LIB_FILE) $(CONFIG_MK) $(SOURCECODE)
	$(MFEM_CXX) $(MFEM_FLAGS) $< -o $@ $(MFEM_LIBS) $(SOURCECODE) $(LDFLAGS)


