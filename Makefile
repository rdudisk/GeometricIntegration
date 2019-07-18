all: reissner

reissner: main.cpp headers/*
	#g++ -o main main.cpp -lgsl -lgslcblas -lm -lpthread -std=gnu++11 
	g++ -o main main.cpp libs/fastgl/fastgl.cpp -lgsl -lgslcblas -lm  -std=c++11

test: test-so3.cpp headers/*
	g++ -o test test-so3.cpp -lgsl -lgslcblas -lm -lpthread -std=c++11

test_dof: test_dof.cpp
	g++ -o test_dof test_dof.cpp 

unit_test: unit_test.cpp
	g++ -o unit_test unit_test.cpp -lgsl -lgslcblas -lm

doxygen:
	doxygen Doxygen

include Makefile.export.NOX

CXX=$(NOX_CXX_COMPILER)
CC=$(NOX_C_COMPILER)
FORT=$(NOX_Fortran_COMPILER)

CXX_FLAGS=$(NOX_CXX_COMPILER_FLAGS)
C_FLAGS=$(NOX_C_COMPILER_FLAGS)
FORT_FLAGS=$(NOX_Fortran_COMPILER_FLAGS)

INCLUDE_DIRS=$(NOX_INCLUDE_DIRS) $(NOX_TPL_INCLUDE_DIRS)
LIBRARY_DIRS=$(NOX_LIBRARY_DIRS) $(NOX_TPL_LIBRARY_DIRS)
LIBRARIES=$(NOX_LIBRARIES) $(NOX_TPL_LIBRARIES)

LINK_FLAGS=$(NOX_EXTRA_LD_FLAGS)

#just assuming that epetra is turned on.
DEFINES=-DMYAPP_EPETRA

midpoint: midpoint.cpp
	$(CXX) $(CXX_FLAGS) -g midpoint.cpp -o midpoint $(LINK_FLAGS) $(INCLUDE_DIRS) $(DEFINES) $(LIBRARY_DIRS) $(LIBRARIES)

midpoint2: midpoint2.cpp
	$(CXX) $(CXX_FLAGS) -g midpoint2.cpp -o midpoint2 $(LINK_FLAGS) $(INCLUDE_DIRS) -I./ $(DEFINES) $(LIBRARY_DIRS) $(LIBRARIES)
