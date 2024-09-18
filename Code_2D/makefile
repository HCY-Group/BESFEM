
CONFIG_MK = $(EBROOTMFEM)/share/mfem/config.mk

include $(CONFIG_MK)

# Remove built-in rule
%: %.cpp

# Replace the default implicit rule for *.cpp files
%: %.cpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	$(MFEM_CXX) $(MFEM_FLAGS) $< -o $@ $(MFEM_LIBS)

docs:
	doxygen ../docs/Doxyfile
