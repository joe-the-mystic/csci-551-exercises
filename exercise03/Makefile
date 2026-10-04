CXX = g++
CXX_FLAGS = -g -Wall -O0 
MPIXX = mpicxx
MPIXX_FLAGS = -g -Wall -O0 


PRODUCT = integrate_seq integrate_par integrate_par_float
CXXFILES = integrate_seq.cpp integrate_par.cpp integrate_par_float.cpp

all: ${PRODUCT}

clean:
	-rm -f *.o *.NEW *~
	-rm -f ${PRODUCT} 

integrate_seq: integrate_seq.cpp
	${CXX} ${CXX_FLAGS} integrate_seq.cpp -o integrate_seq

integrate_par: integrate_par.cpp
	${MPIXX} ${MPIXX_FLAGS} integrate_par.cpp -o integrate_par

integrate_par_float: integrate_par_float.cpp
	${MPIXX} ${MPIXX_FLAGS} integrate_par_float.cpp -o integrate_par_float

