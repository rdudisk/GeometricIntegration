all: testing

doxygen:
	doxygen Doxygen

include /usr/local/trilinos/include/Makefile.export.Epetra

CXX=$(Epetra_CXX_COMPILER)
CC=$(Epetra_C_COMPILER)
FORT=$(Epetra_Fortran_COMPILER)

CXX_FLAGS=$(Epetra_CXX_COMPILER_FLAGS)
C_FLAGS=$(Epetra_C_COMPILER_FLAGS)
FORT_FLAGS=$(Epetra_Fortran_COMPILER_FLAGS)

INCLUDE_DIRS=$(Epetra_INCLUDE_DIRS) $(Epetra_TPL_INCLUDE_DIRS)
LIBRARY_DIRS=$(Epetra_LIBRARY_DIRS) $(Epetra_TPL_LIBRARY_DIRS)
LIBRARIES=$(Epetra_LIBRARIES) $(Epetra_TPL_LIBRARIES)

LINK_FLAGS=$(Epetra_EXTRA_LD_FLAGS)

#just assuming that epetra is turned on.
DEFINES=-DMYAPP_EPETRA

testing: testing.cpp
	$(CXX) $(CXX_FLAGS) -g testing.cpp -o testing $(LINK_FLAGS) $(INCLUDE_DIRS) -I./ $(DEFINES) $(LIBRARY_DIRS) $(LIBRARIES)
