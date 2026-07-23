# General MFEM configuration
#MFEM_DIR ?= ../../..
#MFEM_BUILD_DIR ?= ../../..
#CONFIG_MK = $(MFEM_BUILD_DIR)/config/config.mk
CONFIG_MK = $(MFEM_BUILD_DIR)/share/mfem/config.mk

# Include MFEM's configuration
-include $(CONFIG_MK)

# Compiler and Flags
CXX = mpicxx
CXXFLAGS = -g -O3 -std=c++17 -w

# Use MFEM-provided flags (includes Hypre/Metis/etc)
INCLUDE_FLAGS := $(MFEM_INCFLAGS)
LIB_FLAGS     := $(MFEM_LIBS)

# Source files
SRC_FILES = \
    src/SimulationConfig.cpp \
    src/Initialize_Geometry.cpp \
    src/readtiff.cpp \
    src/Domain_Parameters.cpp \
    src/BoundaryConditions.cpp \
    src/Concentrations_Base.cpp \
    src/ElectrolyteDiffusion.cpp \
    src/ElectrodeCahnHilliard.cpp \
    src/ElectrodeDiffusion.cpp \
    src/Reaction.cpp \
    src/FEMOperators.cpp \
    src/Potentials_Base.cpp \
    src/Adjust.cpp \
    src/Utils.cpp \
    src/dist_solver.cpp \
    src/SimulationState.cpp \
    src/MaterialProperties.cpp \
    src/ElectrodePotential.cpp \
    src/ElectrolytePotential.cpp \
    src/Constants.cpp 

# ====================================

# Output executable
EXEC_DIR = bin
EXEC_NAME = battery_simulation
EXEC = $(EXEC_DIR)/$(EXEC_NAME)
EXEC_SRC_FILES = $(SRC_FILES) src/battery_simulation.cpp


# ================== Test sources ==================
TEST_NAME = unit_tests
TEST_EXEC = $(EXEC_DIR)/$(TEST_NAME)
TEST_SRC_FILES = $(SRC_FILES) tests/unit_tests.cpp tests/unit_tests_main.cpp


# ====================================



# Tiff reading
LDFLAGS = -ltiff #-lCatch2

# Default target
all: $(EXEC)

# Compile the simulation target
$(EXEC): $(EXEC_SRC_FILES)
	@mkdir -p $(EXEC_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) $^ -o $@ $(LIB_FLAGS) $(LDFLAGS)


test: $(TEST_EXEC)
# Build test binary
$(TEST_EXEC): $(TEST_SRC_FILES)
	@mkdir -p $(EXEC_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) $^ -o $@ $(LIB_FLAGS) $(LDFLAGS)


# Clean build artifacts
clean:
	rm -f $(EXEC) *.o *~ *.dSYM *.TVD.*breakpoints
	rm -rf $(OBJ_DIR) $(TEST_BIN)
