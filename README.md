# About
FLPQ is a code to calculate local physical quantities from density matrix obtained by DFT (OpenMX).

# Installation

## makefile
Compiler options (CC, FC, LIB) must be modified for your computer environment.  
If you have already installed OpenMX, you can use the same compiler options.


For example,
```
### ohtaka (CPU:AMD ROME), Intel compiler + Intel MKL + Intel MPI ###
CC = mpiicc -O3 -march=core-avx2 -ip -no-prec-div -qopenmp -I${MKLROOT}/include/fftw -parallel -par-schedule-auto -static-intel -qopenmp-link=static -qopt-malloc-options=3 -qopt-report
FC = mpiifort -O3 -march=core-avx2 -ip -no-prec-div -qopenmp -parallel -par-schedule-auto -static-intel -qopenmp-link=static -qopt-malloc-options=3 -qopt-report
LIB= -L${MKLROOT}/lib/intel64 -qmkl=parallel -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lifcore -lpthread -lm -ldl
```

To compile the code, type
```
$ make
```
Then, executive file flpq is generated.

# Usage
After preparing OpenMX_LPQ.bin and OpenMX_LPQ_DM by the OpenMX calculation,
execute the following command.

```
$ flpq flpq.inp > flpq.log 
```

# Input file example (flpq.inp)
```
dens_origin_unit Ang
dens_origin  -6.0 -6.0 -6.0  # default  1.0e+9 1.0e+9 1.0e+9
Dens_Ngrid  51 51 51         # default  Ngrids are imported from the OpenMX result.
Grid_vectors.Unit   Ang
<Grid_vectors                # Box region for calculating LPQ. default is the lattice vector imported from the OpenMX result.
12.0   0.0   0.0
0.0   12.0   0.0
0.0   0.0   12.0
Grid_vectors>
save_memory_3d   off          # default  off
```


# Contact
Masahiro FUKUDA (ISSP, Univ. of Tokyo)  
masahiro.fukuda__at__issp.u-tokyo.ac.jp  
Please replace `__at__` by @.
