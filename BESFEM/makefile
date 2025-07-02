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
# SRC_FILES = simulation.cpp Constants.cpp \
#             Concentrations_Base.cpp CnP.cpp CnE.cpp \
#             Potentials_Base.cpp PotP.cpp PotE.cpp \
#             Reaction.cpp Current.cpp Initialize_Geometry.cpp \
#             Domain_Parameters.cpp readtiff.cpp

# SRC_FILES = simulation.cpp ../code/Constants.cpp ../code/Initialize_Geometry.cpp readtiff.cpp Domain_Parameters.cpp \
#             Concentrations_Base.cpp CnP.cpp CnE.cpp CnCH.cpp ../code/SolverSteps.cpp  \
#             Potentials_Base.cpp PotP.cpp PotE.cpp Reaction.cpp Current.cpp

# Source files
SRC_FILES = \
    src/simulation.cpp \
    inputs/Constants.cpp \
    src/Initialize_Geometry.cpp \
    src/readtiff.cpp \
    src/Domain_Parameters.cpp \
    src/Concentrations_Base.cpp \
    src/CnE.cpp \
    src/CnCH.cpp \
    src/Reaction.cpp \
    src/SolverSteps.cpp \
    src/Potentials_Base.cpp \
    src/PotP.cpp \
    src/PotE.cpp

# SRC_FILES = Anna_ex_BESFEM_2D_Par_H2_gr.cpp 
            

# Output executable
# EXEC = simulation
# EXEC = serial_simulation
EXEC_DIR = bin
EXEC_NAME = simulation
EXEC = $(EXEC_DIR)/$(EXEC_NAME)

# Tiff reading
LDFLAGS = -ltiff

# Default target
all: $(EXEC)

# Compile the simulation target
$(EXEC): $(SRC_FILES)
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) $^ -o $@ $(LIB_FLAGS) $(LDFLAGS)

# Clean build artifacts
clean:
	rm -f $(EXEC) *.o *~ *.dSYM *.TVD.*breakpoints