# FLPQ

FLPQ calculates local physical quantities from density matrices obtained by
density-functional theory (DFT) calculations using
[OpenMX](https://www.openmx-square.org/). It produces scalar, vector, and
tensor fields, including electronic stress and quantities derived from it.

FLPQ is research software: validate the electronic structure, basis set,
spin configuration, grid spacing, and integration accuracy for the system
of interest.

## Installation

### Prerequisites

- OpenMX output exported as OpenMX_LPQ.bin and OpenMX_LPQ_DM/.
- A C compiler with MPI and OpenMP support.
- Compatible BLAS, LAPACK, and ScaLAPACK libraries, typically Intel MKL.
- Optional HDF5 C library for direct HDF5 output.
- Optional Python 3.9+, NumPy, and h5py for HDF5 conversion and tests.

### Compiler and library configuration

Modify the CC, FC, and LIB options in source/makefile for your computing
environment. **If OpenMX is already installed, its compiler and library
settings are a good starting point**, as in the v0.1.0 documentation.

The following is the Intel compiler + Intel MPI + Intel MKL example used
in the original README for an AMD Rome HPC environment (ohtaka). It is
an example, not a universally portable configuration:

~~~make
CC = mpiicc -O3 -march=core-avx2 -ip -no-prec-div -qopenmp -I${MKLROOT}/include/fftw -parallel -par-schedule-auto -static-intel -qopenmp-link=static -qopt-malloc-options=3 -qopt-report
FC = mpiifort -O3 -march=core-avx2 -ip -no-prec-div -qopenmp -parallel -par-schedule-auto -static-intel -qopenmp-link=static -qopt-malloc-options=3 -qopt-report
LIB = -L${MKLROOT}/lib/intel64 -qmkl=parallel -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lifcore -lpthread -lm -ldl
~~~

Set up your compilers, MPI, and numerical libraries first. Run make from
the source directory; the resulting executable is source/flpq.

Legacy-only build (the default, with no HDF5 dependency):

~~~bash
cd source
make clean
make -j4 HDF5=0
~~~

To enable direct HDF5 output:

~~~bash
cd source
make clean
make -j4 HDF5=1
~~~

Use a clean rebuild whenever switching HDF5=0 and HDF5=1, because the
makefile does not automatically invalidate all object files when this
option changes. Additional HDF5 header/library paths may be needed
on your system.

See [Building](docs/BUILDING.md) for further details.

## Usage

After an OpenMX calculation has generated OpenMX_LPQ.bin and the
OpenMX_LPQ_DM/ directory, prepare a separate working directory:

~~~text
my_calculation/
    OpenMX_LPQ.bin
    OpenMX_LPQ_DM/
    flpq.inp
~~~

FLPQ reads the OpenMX binary files from the **current working directory**.
The input filename is an optional positional argument, defaulting to
flpq.inp.

For example, within compute resources allocated by the scheduler:

~~~bash
mpiexec -n 4 /path/to/FLPQ/source/flpq flpq.inp > flpq.log
~~~

The original v0.1.0 README also illustrated direct execution with one
process:

~~~bash
/path/to/FLPQ/source/flpq flpq.inp > flpq.log
~~~

MPI launcher choice, process binding, Slurm directives and fabric
configuration are site-specific. Do not copy another cluster's
mpiexec or srun settings without checking its MPI integration.

## Input file example (flpq.inp)

The following complete spatial-grid input retains the v0.1.0 example,
including the origin and box vectors:

~~~text
dens_origin.unit Ang
dens_origin  -6.0 -6.0 -6.0
Dens_Ngrid  51 51 51

Grid_vectors.Unit Ang
<Grid_vectors
12.0   0.0   0.0
0.0   12.0   0.0
0.0   0.0   12.0
Grid_vectors>

save_memory_3d off
Output.Format legacy
~~~

Input keyword meanings and defaults:

| Keyword | Meaning | Default |
| --- | --- | --- |
| dens_origin.unit | Unit of the grid origin (e.g. Ang) | See the input parser |
| dens_origin | Grid origin in the specified unit | 1.0e+9 1.0e+9 1.0e+9 (automatic origin) |
| Dens_Ngrid | Number of grid points along the three grid axes | OpenMX-derived grid |
| Grid_vectors.Unit | Unit of the box vectors (e.g. Ang) | See the input parser |
| Grid_vectors block | Three vectors defining the local-quantity box | Lattice vectors from OpenMX |
| save_memory_3d | Enable slice-by-slice, reduced-memory calculation | off |
| Output.Format | Choose legacy, hdf5, or both | legacy |
| Output.Quantities | basic, all, all_scalar, rcp, etc. | all |

**Important correction from v0.1.0:** The previous README used
"dens_origin_unit" (underscores), whereas the actual input reader
recognizes "dens_origin.unit" (with a dot). Use the latter.

If dens_origin, Dens_Ngrid, or the Grid_vectors block is omitted, the
relevant automatic/default behavior applies. The example above
explicitly selects the spatial box and calculation grid.

## Selecting physical quantities

Use `Output.Quantities` in `flpq.inp` to select individual fields or
predefined groups. This works with legacy and HDF5 output:

~~~text
Output.Quantities basic
~~~

Supported presets are `all` (the historical default), `all_scalar`,
`all_vector`, `all_tensor`, `basic` (electron density, spin angular
momentum and noncollinear zeta potential), and `rcp` (kinetic energy
density and chemical potential). The noncollinear-only field in
`basic` is omitted for other spin modes.

See **[the complete list of selectable physical quantities and presets](docs/OUTPUT_QUANTITIES.md)**
for every quantity name, spin restrictions, combination rules, and examples.

## Output formats

### Legacy dat and Cube output

The default Output.Format legacy keeps the original dat/Cube output
behavior and does not need the HDF5 library.

### Optional direct HDF5 output (v0.2.0)

Compile with HDF5=1, then add, for example:

~~~text
Output.Format hdf5
Output.Quantities electron_density,charge_current,stress
Output.HDF5.Precision float64
Output.HDF5.Compression gzip
Output.HDF5.GzipLevel 4
Output.HDF5.OnExisting error
~~~

Output.Format both writes both legacy and HDF5 files. Each HDF5 file
represents one physical quantity and contains its own grid, atom,
units, and component metadata.

For complete keyword choices, compression settings, collision handling,
and limitations, see [HDF5 output](docs/HDF5_OUTPUT.md).

To inspect and convert HDF5 fields to Cube, FieldCube, or legacy dat:

~~~bash
python3 -m pip install -r tools/requirements-hdf5.txt
python3 tools/flpq_convert.py --help
~~~

## Optional symmetric and antisymmetric local stress

Add the following option to the input:

~~~text
Output.Stress.Decomposition on
~~~

The default is off for backward compatibility. When enabled (and the
corresponding quantities are selected), FLPQ writes the full,
symmetric, and antisymmetric local electronic stress tensors:

- LPQ/LPQ.stress.dat
- LPQ/LPQ.stress_S.dat
- LPQ/LPQ.stress_A.dat

The legacy tensor component order is column-major:
xx, yx, zx, xy, yy, zy, xz, yz, zz.

The decomposition is S = (stress + stress^T)/2 and
A = (stress - stress^T)/2. Energy density and local chemical
potential use half the trace of the full stress tensor.

## Regression tests and numerical validation

The methane regression uses a precomputed fixture containing
OpenMX_LPQ.bin, OpenMX_LPQ_DM/, and flpq.inp. It checks output
consistency between 1/2/4 MPI ranks, full-grid versus slice-based
low-memory execution, repeated runs, scalar Cube values, and invalid
or truncated input handling:

~~~bash
python3 tests/regression_methane.py --help
python3 tests/check_stress_decomposition.py --help
~~~

Noncollinear (NC) and spin-orbit coupling (SOC) stress calculations
have also been validated with a two-atom Cr2 system. A previous
fix corrected cross-gradient stress accumulation in the
four-grid-point batched noncollinear paths by using the spin+2
component indices consistently with the scalar paths.

In the corresponding development validation, the Cr2 NC+SOC and
NC-without-SOC cases passed 1/2/4-rank comparisons. The converged
test used an SCF criterion of 1e-7, and grid refinement from
7^3 to 56^3 confirmed pointwise consistency. Integrals involving
spatial derivatives remain sensitive to finite grid resolution.

Historical site-specific Slurm scripts and large OpenMX fixture data
are **not included** in the public release candidate. The portable
test scripts require your own compatible fixtures. Additional HDF5
regressions are described in [Testing](docs/TESTING.md).

## Code organization

- source/main_flpq.c — main control flow and input/output orchestration.
- source/read_input_openmx_lpq.c and source/read_input_DM.c — OpenMX export readers.
- source/Set_Grid.c — local grid and atom-overlap construction.
- source/Set_ddOrbitals_Grid_fuku.c — orbital values and derivatives on grids.
- source/Set_Density_Grid_fuku.c and related modules — density-field evaluation.
- source/Local_physical_quantities.c — physical quantities from grid fields.
- source/out_LPQ.c and source/FLPQ_HDF5.c — legacy/HDF5 output paths.

This release retains the established atom-owned MPI evaluation.
Experimental asynchronous output prefetching, LPT atom redistribution,
and atom-pair/orbital/grid MPI task scheduling are not included because
they did not demonstrate a reliable end-to-end speedup on the tested cases.

## License

[GNU General Public License, version 3](LICENSE).
Follow applicable notices and license terms for any external software.

## Contact

Masahiro FUKUDA (Institute for Solid State Physics,
The University of Tokyo).

For questions and reproducible issue reports, use the
[FLPQ GitHub issue tracker](https://github.com/mfukudaQED/FLPQ/issues).
The previous README included a personal email address; it is
deliberately not reproduced in this public release candidate.
