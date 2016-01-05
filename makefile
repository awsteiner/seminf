FLAGS = -I$(GSL_INC) -I$(O2SCL_INC) -std=c++0x -I$(EIGEN_INC) -I$(HDF5_INC)
LIBS = -L$(GSL_LIB) -L$(O2SCL_LIB) \
	-lo2scl_eos -lo2scl_part -lo2scl_hdf -lo2scl \
	-lreadline -lgsl -lgslcblas -lm

relax.o: relax.cpp relax.h
	g++ $(FLAGS) -c relax.cpp

nr.o: nr.cpp
	g++ $(FLAGS) -c nr.cpp

nr: nr.o relax.o
	g++ $(FLAGS) -o nr nr.o relax.o $(LIBS)

clean:
	rm nr *.o

