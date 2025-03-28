### Compilers r###
CC=g++## for C++
FF=gfortran## for Fortran

###  Directories ###
SRCDIR=${PWD}/src
INCDIR=${PWD}/inc
LHAPDFDIR=$//home/matteo/LHAPDF/include
LHAPDFDIR2=$//home/matteo/LHAPDF/share/LHAPDF
LIBDIR=${PWD}/lib
BLDDIR=${PWD}/bld
GSLDIR=$//home/matteo/gsl-2.8
BOOSTDIR=$//usr/include/boost

### Flags ###
CFLAG=-O3 -Wall -I${INCDIR} -I${LHAPDFDIR} -I${LHAPDFDIR2} -I${GSLDIR} -I${BOOSTDIR}
FFLAG=-O3 -I${INCDIR}

### Paths ###
VPATH=${SRCDIR}

### Source files ###
CFILES=$(wildcard ${SRCDIR}/*.cpp)
FFILES=$(wildcard ${SRCDIR}/*.f)

### Object files ###
COBJS=$(subst .cpp,.o,$(subst ${SRCDIR},${BLDDIR},${CFILES}))
FOBJS=$(subst .f,.o,$(subst ${SRCDIR},${BLDDIR},${FFILES}))

### Libraries ###
LIB=${LIBDIR}/libresum.a
GSLIB=-L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas
STLIB= -L/usr/lib/gcc/x86_64-linux-gnu/13 -lm -lgfortran 
LHAPDFLIB= -L/home/matteo/LHAPDF/lib -lLHAPDF

### Commands ###
all: RUN
lib: ${LIB}

RUN: main.cpp ${LIB}
	${CC} ${CFLAG} -o $@ main.cpp ${LIB} ${STLIB} ${GSLIB} ${LHAPDFLIB} 

${LIB}:	${COBJS} ${FOBJS}
	ar -ruc $@ ${BLDDIR}/*.o

${BLDDIR}%.o: ${SRCDIR}%.f
	cd ${BLDDIR}; ${FF} -c ${FFLAG} $?;

${BLDDIR}%.o: ${SRCDIR}%.cpp
	cd ${BLDDIR}; ${CC} -c ${CFLAG} $<;
	
### Cleaning ###
clean:
	rm -f RUN ${LIBDIR}/*.a ${BLDDIR}/*.o *.log output/*; clear;

