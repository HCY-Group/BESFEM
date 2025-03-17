
# CONFIG_MK = $(EBROOTMFEM)/share/mfem/config.mk

# include $(CONFIG_MK)


# LDFLAGS = -ltiff
# SOURCECODE = ../code/Constants.cpp ../code/Initialize_Geometry.cpp TimeDepOpers.cpp readtiff.cpp MeshMaker.cpp VoxelSolver.cpp VoxelSolver_DG.cpp

# MFEM_CXXFLAGS = -std=c++14 -O3 -DNDEBUG -O2 -ftree-vectorize -march=native -fno-math-errno -fopenmp

# # Remove built-in rule
# %: %.cpp

# # Replace the default implicit rule for *.cpp files
# %: %.cpp $(MFEM_LIB_FILE) $(CONFIG_MK) $(SOURCECODE)
# 	$(MFEM_CXX) $(MFEM_FLAGS) $< -o $@ $(MFEM_LIBS) $(SOURCECODE) $(LDFLAGS)

# General MFEM configuration
MFEM_DIR ?= ../../..
MFEM_BUILD_DIR ?= ../../..
CONFIG_MK = $(MFEM_BUILD_DIR)/config/config.mk

# Include MFEM's configuration
-include $(CONFIG_MK)

# Compiler and Flags
CXX = mpicxx
CXXFLAGS = -g -O3 -std=c++14
INCLUDE_FLAGS = -I$(MFEM_DIR) -I$(MFEM_BUILD_DIR) \
                -I../../../../hypre/src/hypre/include
LIB_FLAGS = -L$(MFEM_BUILD_DIR) -lmfem \
            -L../../../../hypre/src/hypre/lib -lHYPRE \
            -L../../../../metis-4.0 -lmetis -lrt

# Source files for the simulation
SRC_FILES = ../code/Constants.cpp ../code/Initialize_Geometry.cpp TimeDepOpers.cpp readtiff.cpp MeshMaker.cpp VoxelSolver.cpp VoxelSolver_DG.cpp \
			VoxelToDistFunc.cpp

# Output executable
EXEC = voxel

LDFLAGS = -ltiff

# Default target
all: $(EXEC)

# Compile the simulation target
$(EXEC): $(SRC_FILES)
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) $^ -o $@ $(LIB_FLAGS) $(LDFLAGS)

# Clean build artifacts
clean:
	rm -f $(EXEC) *.o *~ *.dSYM *.TVD.*breakpoints

print-vars:
	@echo "MFEM_DIR = $(MFEM_DIR)"
	@echo "MFEM_BUILD_DIR = $(MFEM_BUILD_DIR)"
	@echo "CONFIG_MK = $(CONFIG_MK)"

