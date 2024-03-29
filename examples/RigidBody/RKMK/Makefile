all: RigidBody

include /usr/local/trilinos/include/Makefile.export.NOX

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

RigidBody: RigidBody.cpp
	$(CXX) $(CXX_FLAGS) -g -O3 RigidBody.cpp -o RigidBody $(LINK_FLAGS) $(INCLUDE_DIRS) -I./../../../ $(DEFINES) $(LIBRARY_DIRS) $(LIBRARIES)

figure:
	python ../../../visualisation/rigid-body/projection.py results.csv
	cp -t . ../../../visualisation/rigid-body/exact0.csv ../../../visualisation/rigid-body/exact1.csv
	sed 's/CSV_FILENAME/output.csv/g' ../../../visualisation/rigid-body/template.tex > tmp.tex
	latexmk -f -pdf -pdflatex="pdflatex -shell-escape -interaction=nonstopmode" -use-make tmp.tex
	pdf2svg tmp.pdf figure.svg
	gs -q -dNOCACHE -dNOPAUSE -dBATCH -dSAFER -sDEVICE=eps2write -sOutputFile=figure.eps tmp.pdf
	rm output.csv exact0.csv exact1.csv tmp.*
