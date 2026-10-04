CXX = g++
CXX_FLAGS = -g -Wall -O0 
MPIXX = mpicxx
MPIXX_FLAGS = -g -Wall -O0 


PRODUCT = integrate_seq integrate_par
CXXFILES = integrate_seq.cpp integrate_par.cpp

all: ${PRODUCT}

clean:
	-rm -f *.o *.NEW *~
	-rm -f ${PRODUCT} 

integrate_seq: integrate_seq.cpp
	${CXX} ${CXX_FLAGS} integrate_seq.cpp -o integrate_seq

integrate_par: integrate_par.cpp
	${MPIXX} ${MPIXX_FLAGS} integrate_par.cpp -o integrate_par
