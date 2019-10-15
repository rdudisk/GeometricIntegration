all: testing

doxygen:
	doxygen Doxygen

include /usr/local/trilinos/include/Makefile.export.Trilinos
include /usr/local/trilinos/include/Makefile.export.NOX
include /usr/local/trilinos/include/Makefile.export.Epetra

CXX=$(Trilinos_CXX_COMPILER)
CC=$(Trilinos_C_COMPILER)
#FORT=$(Epetra_Fortran_COMPILER)

CXX_FLAGS=$(Trilinos_CXX_COMPILER_FLAGS)
C_FLAGS=$(Trilinos_C_COMPILER_FLAGS)
#FORT_FLAGS=$(Epetra_Fortran_COMPILER_FLAGS)

INCLUDE_DIRS=$(Epetra_INCLUDE_DIRS) $(Epetra_TPL_INCLUDE_DIRS) $(NOX_INCLUDE_DIRS) $(NOX_TPL_INCLUDE_DIRS)
LIBRARY_DIRS=$(Epetra_LIBRARY_DIRS) $(Epetra_TPL_LIBRARY_DIRS) $(NOX_LIBRARY_DIRS) $(NOX_TPL_LIBRARY_DIRS)
LIBRARIES=$(Epetra_LIBRARIES) $(Epetra_TPL_LIBRARIES) $(NOX_LIBRARIES) $(NOX_TPL_LIBRARIES)

LINK_FLAGS=$(Epetra_EXTRA_LD_FLAGS) $(NOX_EXTRA_LD_FLAGS)

#just assuming that epetra and NOX are turned on.
DEFINES=-DMYAPP_EPETRA -DMYAPP_NOX

testing: testing.cpp
	$(CXX) $(CXX_FLAGS) -g testing.cpp -o testing $(LINK_FLAGS) $(INCLUDE_DIRS) -I./ $(DEFINES) $(LIBRARY_DIRS) $(LIBRARIES)
