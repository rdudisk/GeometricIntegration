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

solveur: tmp/solveur.o tmp/myproblem.o tmp/legendre.o tmp/io.o tmp/pugixml.o
	$(CXX) $(CXX_FLAGS) -g tmp/solveur.o tmp/myproblem.o tmp/legendre.o tmp/io.o tmp/pugixml.o -o solveur $(LINK_FLAGS) $(LIBRARY_DIRS) $(LIBRARIES) -lgsl -lgslcblas -lm
tmp/solveur.o: solveur.cpp header.hpp
	$(CXX) $(CXX_FLAGS) -c solveur.cpp -o tmp/solveur.o $(LINK_FLAGS) $(INCLUDE_DIRS) -I/home/ttr2/these/ressources/code/geomi $(DEFINES)
tmp/myproblem.o: myproblem.cpp header.hpp
	$(CXX) $(CXX_FLAGS) -c myproblem.cpp -o tmp/myproblem.o $(LINK_FLAGS) $(INCLUDE_DIRS) -I/home/ttr2/these/ressources/code/geomi $(DEFINES)
tmp/legendre.o: legendre.cpp header.hpp
	$(CXX) $(CXX_FLAGS) -c legendre.cpp -o tmp/legendre.o $(LINK_FLAGS) $(INCLUDE_DIRS) -I/home/ttr2/these/ressources/code/geomi $(DEFINES)
tmp/io.o: io.cpp header.hpp
	$(CXX) $(CXX_FLAGS) -c io.cpp -o tmp/io.o $(LINK_FLAGS) $(INCLUDE_DIRS) -I/home/ttr2/these/ressources/code/geomi $(DEFINES)
tmp/pugixml.o: pugixml.cpp
	g++ -O3 -c pugixml.cpp -o tmp/pugixml.o

test: solveur
	date && ./solveur > log.txt && date

animation:
	echo "Creating plots... "
	python anim/plot-reissner.py results.csv anim/img
	echo "done\n"
	echo "Creating movie..."
	ffmpeg -y -framerate 25 -i anim/img/imgReissner%d.png movie.mp4
	echo "done"

clean:
	/bin/rm -f tmp/solveur.o tmp/classes.o

mrproper: clean
	/bin/rm solveur log.txt results.csv movie.avi
