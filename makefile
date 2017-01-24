ifeq ($(USER),awsteiner)

else

GSL_INC = .
O2SCL_INC = .
EIGEN_INC = .
HDF5_INC = .
GSL_LIB = .
O2SCL_LIB = .
HDF5_LIB = .

endif

FLAGS = -I$(GSL_INC) -I$(O2SCL_INC) -std=c++0x -I$(EIGEN_INC) \
	-I$(HDF5_INC) -Wno-deprecated-declarations
#-DUSE_EIGEN

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

check:
	nr check > check_nr.scr
	tail -n 2 check_nr.scr
	rel check > check_rel.scr
	tail -n 2 check_rel.scr

doc: empty
	git rev-parse HEAD | awk \
		'{print "<a href=\"http://github.com/awsteiner/seminf/tree/" $$1 "\">" $$1 "</a>"}' \
		 > doc/rev.txt
	cd doc; doxygen doxyfile

empty:
