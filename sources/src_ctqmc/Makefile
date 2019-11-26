include ./Makefile.in

#C++ = mpicxx_gnu
#CFLAGS = -O3 -D_MPI -I/opt/mpich-1.2.6/include #-D_DEBUG #-DNO_ARG_CHECK
#LIBS = -static -L/opt/pgi/linux86/6.0/lib -llapack -lblas -lpgftnrtl -lpgc -lgsl -lgslcblas

#PC++ = mpicxx # g++
#CFLAGS =  -O3 -D_MPI -I/opt/mpich-1.2.7-gnu/include #-DAS # -g -D_DEBUG -Wall 
#LIBS = -L/opt/intel/mkl/10.0.1.014/lib/32 -lmkl_lapack -lmkl -lguide -lgsl -pthread
#LIBS = -i-static -L/opt/pgi/linux86-64/7.0-3/lib -llapack -lblas -lpgftnrtl -lpgc -lgsl -lgslcblas -lrt

CFLAGS = $(PFLAGS)  -DAS #-Ofast -xAVX -fpermissive # -D_TIME #-D_LOGGING #-DAS 
LIBS = $(PLIBS)

GHEADERS = assert.h complex.h random.h sblas.h sfunction.h smesh.h sutil.h zeroin.h
QHEADERS = common.h  inout.h intervals.h local.h matrixm.h mpi.h bcast.h number.h operators.h state.h stateim.h segment.h

all : ctqmc 

ctqmc : ctqmc.o SMatrix1.o 
	$(PC++) $(CFLAGS) -o $@ ctqmc.o SMatrix1.o $(LIBS)

SMatrix1.o : SMatrix1.cc sfunction.h
	$(PC++) -c $(CFLAGS) SMatrix1.cc 

ctqmc.o : ctqmc.cc $(GHEADERS) $(QHEADERS) 
	$(PC++) $(CFLAGS) -c ctqmc.cc

clean :
	rm -f ctqmc.o ctqmc SMatrix1.o

.SUFFIXES : .cc
.cc.o:
	$(PC++) $(CFLAGS) -c $<

.SUFFIXES : .f
.f.o:
	$(F77) $(FFLAGS) -c $<





