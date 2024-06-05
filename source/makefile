
### ohtaka, Intel compiler + Intel MKL + Open MPI
### used module ###
#module load oneapi_compiler/2023.0.0 oneapi_mkl/2023.0.0 openmpi/4.1.5-oneapi-2023.0.0-classic
#CC = mpicc -O3 -march=core-avx2 -ip -no-prec-div -qopenmp -I${MKLROOT}/include/fftw -parallel -par-schedule-auto -static-intel -qopenmp-link=static -qopt-malloc-options=3 -qopt-report
#FC = mpifort -O3 -march=core-avx2 -ip -no-prec-div -qopenmp -parallel -par-schedule-auto -static-intel -qopenmp-link=static -qopt-malloc-options=3 -qopt-report
#LIB= -L${MKLROOT}/lib/intel64 -mkl=parallel -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi  -lifcoremt

### ohtaka, Intel compiler + Intel MKL + Intel MPI
### used module ###
#module load oneapi_compiler/2023.0.0 oneapi_mkl/2023.0.0 oneapi_mpi/2023.0.0
CC = mpiicc -O3 -march=core-avx2 -ip -no-prec-div -qopenmp -I${MKLROOT}/include/fftw -parallel -par-schedule-auto -static-intel -qopenmp-link=static -qopt-malloc-options=3 -qopt-report
FC = mpiifort -O3 -march=core-avx2 -ip -no-prec-div -qopenmp -parallel -par-schedule-auto -static-intel -qopenmp-link=static -qopt-malloc-options=3 -qopt-report
LIB= -L${MKLROOT}/lib/intel64 -qmkl=parallel -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lifcore -lpthread -lm -ldl

#CC = mpicc -O3 -xHOST -qopenmp -ip -no-prec-div -mkl=parallel -I$(MKLROOT)/include/fftw -shared-intel -pg
#FC = mpif90 -O3 -xHOST -qopenmp -ip -no-prec-div -mkl=parallel -I$(MKLROOT)/include/fftw -shared-intel -pg
#LIB = -L$(MKLROOT)/lib/intel64 -lmpi_mpifh -lifcore

##
## macloud
##
#CC = mpicc -O3 -xHOST -ip -Dscalapack -no-prec-div -qopenmp -I${MKLROOT}/include/fftw
#FC = mpif90 -O3 -xHOST -ip -Dscalapack -no-prec-div -qopenmp
#LIB= -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_openmpi_lp64 -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -liomp5 -lpthread -lm -ldl
#
## 
## Oakforest-PACS
## 
#CC = mpiicc   -O3 -axMIC-AVX512 -qopenmp -ip -no-prec-div -mkl=parallel -I$(MKLROOT)/include/fftw -shared-intel
#FC = mpiifort -O3 -axMIC-AVX512 -qopenmp -ip -no-prec-div -mkl=parallel -I$(MKLROOT)/include/fftw -shared-intel
#LIB = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lifcore -lmpi

CFLAGS  = -g 

OBJS    = main_flpq.o \
          Input_std.o \
		  Inputtools.o \
          read_input_openmx_lpq.o\
          read_input_DM.o\
          util.o\
          Set_Grid.o\
          Get_Cnt_ddOrbitals.o\
          Get_ddOrbitals.o\
          Set_ddOrbitals_Grid_fuku.o\
          Set_Density_Grid_fuku.o\
          Set_dDensity_Grid_fuku.o\
          Set_ddDensity_Grid_fuku.o\
          Set_dd2Density_Grid_fuku.o\
          Set_Density_Grid_fuku_i.o\
          Set_dDensity_Grid_fuku_i.o\
          Set_ddDensity_Grid_fuku_i.o\
          Set_dd2Density_Grid_fuku_i.o\
          out_LPQ.o\
          Local_physical_quantities.o\
		  Set_Allocate_Atom2CPU.o\
          Set_Inf_SndRcv.o\
		  lapack_dstegr1.o\
          Set_Grid_Origin_N.o\
		  Free_Arrays.o\
		  Free_Arrays_in_Set_Grid.o\
          #lapack_prototype.o\
          #xyz2spherical.o\
          #Get_ddOrbitals_QuinticSpline.o\
          #read_input_flpq.o\

# PROG    = openmx.exe
# PROG    = openmx

#
# set program name
# destination directory
#

PROG    = flpq
DESTDIR = ../work

# flpq

flpq:	$(OBJS)
	$(CC) $(OBJS) $(STACK) $(LIB) -lm -o flpq

main_flpq.o: main_flpq.c flpq.h
	$(CC) -c main_flpq.c
Input_std.o: Input_std.c flpq.h Inputtools.h
	$(CC) -c Input_std.c
Inputtools.o: Inputtools.c flpq.h Inputtools.h
	$(CC) -c Inputtools.c
#read_input_flpq.o: read_input_flpq.c flpq.h
#	$(CC) -c read_input_flpq.c
read_input_openmx_lpq.o: read_input_openmx_lpq.c flpq.h
	$(CC) -c read_input_openmx_lpq.c
read_input_DM.o: read_input_DM.c flpq.h
	$(CC) -c read_input_DM.c
util.o: util.c flpq.h
	$(CC) -c util.c
Set_Grid.o: Set_Grid.c flpq.h
	$(CC) -c Set_Grid.c
Get_Cnt_ddOrbitals.o: Get_Cnt_ddOrbitals.c flpq.h
	$(CC) -c Get_Cnt_ddOrbitals.c
Get_ddOrbitals.o: Get_ddOrbitals.c flpq.h
	$(CC) -c Get_ddOrbitals.c
Set_ddOrbitals_Grid_fuku.o: Set_ddOrbitals_Grid_fuku.c flpq.h
	$(CC) -c Set_ddOrbitals_Grid_fuku.c
Set_Density_Grid_fuku.o: Set_Density_Grid_fuku.c flpq.h
	$(CC) -c Set_Density_Grid_fuku.c
Set_dDensity_Grid_fuku.o: Set_dDensity_Grid_fuku.c flpq.h
	$(CC) -c Set_dDensity_Grid_fuku.c
Set_ddDensity_Grid_fuku.o: Set_ddDensity_Grid_fuku.c flpq.h
	$(CC) -c Set_ddDensity_Grid_fuku.c
Set_dd2Density_Grid_fuku.o: Set_dd2Density_Grid_fuku.c flpq.h
	$(CC) -c Set_dd2Density_Grid_fuku.c
Set_Density_Grid_fuku_i.o: Set_Density_Grid_fuku_i.c flpq.h
	$(CC) -c Set_Density_Grid_fuku_i.c
Set_dDensity_Grid_fuku_i.o: Set_dDensity_Grid_fuku_i.c flpq.h
	$(CC) -c Set_dDensity_Grid_fuku_i.c
Set_ddDensity_Grid_fuku_i.o: Set_ddDensity_Grid_fuku_i.c flpq.h
	$(CC) -c Set_ddDensity_Grid_fuku_i.c
Set_dd2Density_Grid_fuku_i.o: Set_dd2Density_Grid_fuku_i.c flpq.h
	$(CC) -c Set_dd2Density_Grid_fuku_i.c
out_LPQ.o: out_LPQ.c flpq.h Local_physical_quantities.h
	$(CC) -c out_LPQ.c
Local_physical_quantities.o: Local_physical_quantities.c flpq.h Local_physical_quantities.h
	$(CC) -c Local_physical_quantities.c
Set_Allocate_Atom2CPU.o: Set_Allocate_Atom2CPU.c flpq.h
	$(CC) -c Set_Allocate_Atom2CPU.c
Set_Inf_SndRcv.o: Set_Inf_SndRcv.c flpq.h
	$(CC) -c Set_Inf_SndRcv.c
lapack_dstegr1.o: lapack_dstegr1.c flpq.h
	$(CC) -c lapack_dstegr1.c
Set_Grid_Origin_N.o: Set_Grid_Origin_N.c flpq.h
	$(CC) -c Set_Grid_Origin_N.c
Free_Arrays.o: Free_Arrays.c flpq.h
	$(CC) -c Free_Arrays.c
Free_Arrays_in_Set_Grid.o: Free_Arrays_in_Set_Grid.c flpq.h
	$(CC) -c Free_Arrays_in_Set_Grid.c

#.o: .c flpq.h
#	$(CC) -c .c

clean:
	rm -f $(PROG) $(OBJS) *.optrpt
