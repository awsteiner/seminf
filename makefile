# ----------------------------------------------------------------------
# This section has variables which may need to be modified. The rest
# of the makefile should not need too much editing.
# ----------------------------------------------------------------------

# This variable should include the directories for the O2scl, GSL, and
# HDF5 libraries. In my configuration, I use the environment variables
# O2SCL_LIB, but you can just replace this entire line with whatever
# you need.
LIB_DIRS = $(LDFLAGS)

# -L$(O2SCL_LIB) -L$(HDF5_LIB)

# This variable may need to be modified to specify the include
# directories for the GSL, boost, HDF5, and O2scl header files. If
# O2scl was installed with Eigen or armadillo support, those header
# directories may need to be here also.
INC_DIRS = $(CXXFLAGS)

# -I$(O2SCL_INC) -I$(HDF5_INC)

# Generic (no MPI necessary) C++ compiler (e.g. g++)
# CXX = 

# Set this variables to empty if you do not have GNU readline support
READLINE_LIBS = -lreadline -lncurses

# Basic compiler flags
COMPILER_FLAGS = -std=c++0x -O3

# ----------------------------------------------------------------------
# Secondary variables:
# ----------------------------------------------------------------------

ALL_FLAGS = $(COMPILER_FLAGS) $(INC_DIRS) $(READLINE_VAR) 

LIBS = -lo2scl_hdf -lo2scl_eos -lo2scl_part -lo2scl \
	-lhdf5_hl -lhdf5 -lgsl -lgslcblas -lm $(READLINE_LIBS)

# ----------------------------------------------------------------------
# Targets:
# ----------------------------------------------------------------------

help:
	@echo "nr: Non-relativistic version "
	@echo "rel: Relativsitic version"
	@echo "clean: "
	@echo "check: "
	@echo "doc: "
	@echo "sync-doc: "

nr.o: nr.cpp
	$(CXX) $(ALL_FLAGS) -c nr.cpp

nr: nr.o
	$(CXX) $(ALL_FLAGS) -o nr nr.o $(LIB_DIRS) $(LIBS)

rel.o: rel.cpp
	$(CXX) $(ALL_FLAGS) -c rel.cpp

rel: rel.o
	$(CXX) $(ALL_FLAGS) -o rel rel.o $(LIB_DIRS) $(LIBS)

clean:
	rm -f nr rel *.o

check:
	./nr -check > check_nr.scr
	tail -n 2 check_nr.scr
	./rel -check > check_rel.scr
	tail -n 2 check_rel.scr

# This version shows output to the screen which ensures travis doesn't
# time out
travis-check:
	./nr -check 
	./rel -check

doc: empty
# Copy most recent tag files
	cd doc; cp ~/o2scl/doc/o2scl/o2scl.tag .
	cd doc; cp ~/o2scl/doc/o2scl/part/o2scl_part.tag .
	cd doc; cp ~/o2scl/doc/o2scl/eos/o2scl_eos.tag .
# Get most recent commit hash
	git rev-parse HEAD | awk \
		'{print "`" $$1 " <http://github.com/awsteiner/seminf/tree/" $$1 ">`_"}' \
		 > sphinx/commit.rst
# Parse bibliography
	cd sphinx/static; cat bib_header.txt > ../bib.rst
	cd sphinx/static; btmanip -parse seminf.bib -rst ../bib_temp.rst
	cd sphinx; cat bib_temp.rst >> bib.rst; rm -f bib_temp.rst
# Run Doxygen
	cd doc; doxygen doxyfile
# Run sphinx
	cd sphinx; make html

sync-doc:
	rsync -Cavzu sphinx/build/html/* $(STATIC_DOC_DIR)/seminf

empty:
