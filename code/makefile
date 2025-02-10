# General MFEM configuration
MFEM_DIR ?= ../..
MFEM_BUILD_DIR ?= ../..
CONFIG_MK = $(MFEM_BUILD_DIR)/config/config.mk

# Include MFEM's configuration
-include $(CONFIG_MK)

# Compiler and Flags
CXX = mpicxx
CXXFLAGS = -g -O3 -std=c++14
INCLUDE_FLAGS = -I$(MFEM_DIR) -I$(MFEM_BUILD_DIR) \
                -I../../../hypre/src/hypre/include
LIB_FLAGS = -L$(MFEM_BUILD_DIR) -lmfem \
            -L../../../hypre/src/hypre/lib -lHYPRE \
            -L../../../metis-4.0 -lmetis -lrt

# Source files for the simulation
SRC_FILES = simulation.cpp Constants.cpp Initialize_Geometry.cpp \
            readtiff.cpp

# Output executable
EXEC = simulation

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