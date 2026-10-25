
CXX = g++
CXX_FLAGS = -g -Wall -O0 
MPIXX = mpicxx
MPIXX_FLAGS = -g -Wall -O0 


PRODUCT = interpolate integrate default_profile.o default_profile.h
CXXFILES = interpolate.cpp integrate.cpp
HFILE = accel_1hz.h

all: ${PRODUCT}

clean:
	-rm -f *.o *.NEW *~
	-rm -f ${PRODUCT} 

interpolate: interpolate.cpp
	${CXX} ${CXX_FLAGS} interpolate.cpp -o interpolate

integrate: integrate.cpp default_profile.o
	${MPIXX} ${MPIXX_FLAGS} integrate.cpp -o integrate

default_profile.o: ${HFILE}
	cp ${HFILE} default_profile.h
	${CXX} ${CXX_FLAGS} default_profile.h -o default_profile.o
