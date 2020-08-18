all: solveur

include /usr/local/trilinos/include/Makefile.export.Trilinos
include /usr/local/trilinos/include/Makefile.export.NOX

CXX=$(Trilinos_CXX_COMPILER)
CC=$(Trilinos_C_COMPILER)

CXX_FLAGS=$(Trilinos_CXX_COMPILER_FLAGS)
C_FLAGS=$(Trilinos_C_COMPILER_FLAGS)

INCLUDE_DIRS=$(NOX_INCLUDE_DIRS) $(NOX_TPL_INCLUDE_DIRS)
LIBRARY_DIRS=$(NOX_LIBRARY_DIRS) $(NOX_TPL_LIBRARY_DIRS)
LIBRARIES=$(NOX_LIBRARIES) $(NOX_TPL_LIBRARIES)

LINK_FLAGS=$(NOX_EXTRA_LD_FLAGS)

DEFINES= -DMYAPP_NOX

solveur: solveur.o classes.o
	$(CXX) $(CXX_FLAGS) -g solveur.o classes.o -o solveur $(LINK_FLAGS) $(LIBRARY_DIRS) $(LIBRARIES) -lgsl -lgslcblas -lm
solveur.o: solveur.cpp classes.hpp
	$(CXX) $(CXX_FLAGS) -c solveur.cpp -o solveur.o $(LINK_FLAGS) $(INCLUDE_DIRS) -I/home/ttr2/these/ressources/code/reissner-reborn/Geomi $(DEFINES)
classes.o: classes.cpp classes.hpp
	$(CXX) $(CXX_FLAGS) -c classes.cpp -o classes.o $(LINK_FLAGS) $(INCLUDE_DIRS) -I/home/ttr2/these/ressources/code/reissner-reborn/Geomi $(DEFINES)

test: solveur
	date && ./solveur > log.txt && date

clean:
	/bin/rm *.o
