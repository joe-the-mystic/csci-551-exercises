INCLUDE_DIRS = 
LIB_DIRS = 
CC = gcc
CXX = g++

CDEFS=
CFLAGS= -O2 -g -Wall -fopenmp $(INCLUDE_DIRS) $(CDEFS)
CFLAGS2= -O0 -g -Wall -fopenmp $(INCLUDE_DIRS) $(CDEFS)
LIBS=

PRODUCT= sumnums_omp_simple sumnums_omp find_primes copy_pgm rotate_pgm rotate_pgms

HFILES= 
CFILES= sumnums_omp_simple.c sumnums_omp.c find_primes.cpp copy_pgm.cpp rotate_pgm.cpp rotate_pgms.cpp

SRCS= ${HFILES} ${CFILES}
OBJS= ${CFILES:.c=.o}

all:	${PRODUCT}

clean:
	-rm -f *.o *.NEW *~
	-rm -f ${PRODUCT} ${DERIVED} ${GARBAGE}

sumnums_omp_simple:	sumnums_omp_simple.c
	$(CC) $(CFLAGS2) -fopenmp -o $@ sumnums_omp_simple.c 

sumnums_omp:	sumnums_omp.c
	$(CC) $(CFLAGS2) -fopenmp -o $@ sumnums_omp.c

find_primes:	find_primes.cpp
	$(CXX) $(CFLAGS2) -fopenmp -o $@ find_primes.cpp

copy_pgm:	copy_pgm.cpp
	$(CXX) $(CFLAGS2) -fopenmp -o $@ copy_pgm.cpp

rotate_pgm:	rotate_pgm.cpp
	$(CXX) $(CFLAGS2) -fopenmp -o $@ rotate_pgm.cpp

rotate_pgms:	rotate_pgms.cpp
	$(CXX) $(CFLAGS2) -fopenmp -o $@ rotate_pgms.cpp
