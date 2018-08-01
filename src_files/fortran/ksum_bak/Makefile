
F77 = mpiifort
PFLAGS = -D_MPI -DMPICH_IGNORE_CXX_SEEK
CFLAGS = $(PFLAGS) -O3 -DAS -msse3 #-xSSE4.2
FFLAGS = -O2 #-xhost

#MKLROOT=/export/share/compilers/intel/composer_xe_2013_sp1.1.106/mkl
PLIBS=-L$(MKLROOT)/lib/intel64/ \
    -Wl,--start-group \
    -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
    -Wl,--end-group -lpthread -lm
LIBS = $(PLIBS)

#all : cmp_mu cmp_mu_sp dmft_ksum dmft_ksum_up dmft_ksum_dn
all : dmft_ksum_sp dmft_ksum 

cmp_mu : cmp_mu.f
	$(F77) $(FFLAGS) -o $@ cmp_mu.f $(LIBS)

dmft_ksum_sp : dmft_ksum_sp.f
	$(F77) $(FFLAGS) -o $@ dmft_ksum_sp.f $(LIBS)
dmft_ksum : dmft_ksum.f
	$(F77) $(FFLAGS) -o $@ dmft_ksum.f $(LIBS)

clean :
	rm -f cmp_mu dmft_ksum dmft_ksum_sp 


.SUFFIXES : .cc
.cc.o:
	$(PC++) $(CFLAGS) -c $<

.SUFFIXES : .f
.f.o:
	$(F77) $(FFLAGS) -c $<

