# Building FLPQ

FLPQ is a hybrid MPI/OpenMP C application with numerical-library
dependencies. The supplied source/makefile defaults to an Intel compiler,
Intel MPI, and Intel MKL configuration; it is an **example**, not an
automatically portable configuration.

## Requirements

1. MPI C compiler and compatible MPI runtime.
2. OpenMP-capable C compiler.
3. BLAS, LAPACK and ScaLAPACK libraries compatible with the MPI implementation.
4. Optional HDF5 C headers and library for HDF5=1.

A working OpenMX installation can often provide a useful starting point for
the MPI and numerical library flags, but the OpenMX build does not need to be
modified to compile FLPQ.

## Minimal build

Adjust source/makefile for your system, then:

~~~bash
cd source
make clean
make -j4 HDF5=0
~~~

The executable is source/flpq. Only HDF5=1 links the HDF5 writer; HDF5=0
links the small stub and continues to use the established legacy output.

## HDF5-enabled build

Make sure HDF5 headers and libraries can be found by the C compiler/linker.
The current makefile appends -lhdf5 when HDF5=1, and may need environment-
specific HDF5 include and library flags.

~~~bash
cd source
make clean
make -j4 HDF5=1
~~~

Important: make clean is necessary whenever switching HDF5=0 and HDF5=1.
The existing makefile does not track the selected writer as a dependency of
the executable; stale objects may otherwise be reused. Do not run clean in
a directory that contains uncommitted binaries or objects you need to keep.

## Execution

Execute in a calculation directory with OpenMX_LPQ.bin,
OpenMX_LPQ_DM/, and flpq.inp. Reserve computing resources first; the
MPI launcher, module commands and fabric settings are system-specific.

~~~bash
mpiexec -n 4 /absolute/path/to/FLPQ/source/flpq flpq.inp
~~~

The optional positional input filename defaults to flpq.inp. The OpenMX
binary files must be in the current calculation directory.

## Supported scope

This release candidate was validated previously with Intel MPI/MKL on an
HPC cluster using precomputed methane and Cr2 fixtures. This does **not**
establish portability to all GCC/OpenMPI or vendor BLAS configurations.
A full compile and regression test is required on each new build platform.
