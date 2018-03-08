
help:
	@echo "nr: "
	@echo "rel: "
	@echo "clean: "
	@echo "check: "
	@echo "doc: "
	@echo "sync-doc: "

ifeq ($(USER),awsteiner)

else

GSL_INC = .
O2SCL_INC = /usr/local/include
EIGEN_INC = .
HDF5_INC = .
GSL_LIB = .
O2SCL_LIB = /usr/local/lib
HDF5_LIB = .

endif

FLAGS = -I$(GSL_INC) -I$(O2SCL_INC) -std=c++0x -I$(EIGEN_INC) \
	-I$(HDF5_INC) -Wno-deprecated-declarations
#-DUSE_EIGEN

LIBS = -L$(GSL_LIB) -L$(O2SCL_LIB) -L$(HDF5_LIB) \
	-lo2scl_hdf -lo2scl_eos -lo2scl_part -lo2scl -lhdf5 -lhdf5_hl \
	-lreadline -lgsl -lgslcblas -lm

nr.o: nr.cpp
	$(CXX) $(FLAGS) -c nr.cpp

nr: nr.o
	$(CXX) $(FLAGS) -o nr nr.o $(LIBS)

rel.o: rel.cpp
	$(CXX) $(FLAGS) -c rel.cpp

rel: rel.o
	$(CXX) $(FLAGS) -o rel rel.o $(LIBS)

clean:
	rm -f nr *.o rel

check:
	./nr check > check_nr.scr
	tail -n 2 check_nr.scr
	./rel check > check_rel.scr
	tail -n 2 check_rel.scr

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
	sudo cp -r sphinx/build/html/* $(STATIC_DOC_DIR)/seminf

empty:
