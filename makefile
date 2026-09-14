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
    src/battery_simulation.cpp \
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

EXEC_DIR = bin
EXEC_NAME = battery_simulation
EXEC = $(EXEC_DIR)/$(EXEC_NAME)


# ================== Test sources ==================
TEST_DIR    := tests
TEST_SRC    := $(TEST_DIR)/besfem_oop_tests.cpp
TEST_BIN    := $(EXEC_DIR)/besfem_oop_tests

# Reuse your sources but drop simulation.cpp (it has its own main()).
TEST_SRC_FILES := $(filter-out src/battery_simulation.cpp,$(SRC_FILES))

OBJ_DIR   := build
TEST_OBJS := $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(TEST_SRC_FILES))

# ====================================

# Tiff reading
LDFLAGS = -ltiff

# Default target
all: $(EXEC)

# Compile the simulation target
$(EXEC): $(SRC_FILES)
	@mkdir -p $(EXEC_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) $^ -o $@ $(LIB_FLAGS) $(LDFLAGS)


# Pattern rule for test objects
$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -c $< -o $@

# Build test binary
$(TEST_BIN): $(TEST_OBJS) $(TEST_SRC)
	@mkdir -p $(EXEC_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) $(TEST_OBJS) $(TEST_SRC) -o $@ $(LIB_FLAGS) $(LDFLAGS)


# Convenience targets
tests: $(TEST_BIN)
	@mkdir -p $(EXEC_DIR)
# 	mpirun -np 1 $(TEST_BIN)

test: tests

# Clean build artifacts
clean:
	rm -f $(EXEC) *.o *~ *.dSYM *.TVD.*breakpoints
	rm -rf $(OBJ_DIR) $(TEST_BIN)
