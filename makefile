FLAGS = -I$(GSL_INC) -I$(O2SCL_INC) -std=c++0x -I$(EIGEN_INC) -I$(HDF5_INC) \
	-Wno-deprecated-declarations

LIBS = -L$(GSL_LIB) -L$(O2SCL_LIB) -L$(HDF5_LIB) \
	-lo2scl_eos -lo2scl_part -lo2scl_hdf -lo2scl -lhdf5 -lhdf5_hl \
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

doc: empty
	git rev-parse HEAD | awk \
		'{print "<a href=\"http://github.com/awsteiner/seminf/tree/" $$1 "\">" $$1 "</a>"}' \
		 > doc/rev.txt
	cd doc; doxygen doxyfile

empty:
