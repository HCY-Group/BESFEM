
# Use the MFEM build directory
MFEM_DIR ?= ../../mfem-4.7/
MFEM_BUILD_DIR ?= ../../mfem-4.7/
CONFIG_MK = $(MFEM_BUILD_DIR)/config/config.mk
# Use the MFEM install directory
# MFEM_INSTALL_DIR = ../mfem
# CONFIG_MK = $(MFEM_INSTALL_DIR)/share/mfem/config.mk

MFEM_LIB_FILE = mfem_is_not_built
-include $(CONFIG_MK)


LDFLAGS = -ltiff
SOURCECODE = TimeDepOpers.cpp readtiff.cpp


# Remove built-in rule
%: %.cpp

# Replace the default implicit rule for *.cpp files
%: %.cpp $(MFEM_LIB_FILE) $(CONFIG_MK) $(SOURCECODE)
	$(MFEM_CXX) $(MFEM_FLAGS) $< -o $@ $(MFEM_LIBS) $(SOURCECODE) $(LDFLAGS)

all: BESFEM_2D_Par VoxelToDistFunc
