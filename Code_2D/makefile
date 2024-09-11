
CONFIG_MK = $(EBROOTMFEM)/share/mfem/config.mk

include $(CONFIG_MK)


LDFLAGS = -ltiff
SOURCECODE = TimeDepOpers.cpp readtiff.cpp MeshMaker.cpp VoxelSolver.cpp



# Remove built-in rule
%: %.cpp

# Replace the default implicit rule for *.cpp files
%: %.cpp $(MFEM_LIB_FILE) $(CONFIG_MK) $(SOURCECODE)
	$(MFEM_CXX) $(MFEM_FLAGS) $< -o $@ $(MFEM_LIBS) $(SOURCECODE) $(LDFLAGS)


