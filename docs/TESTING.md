# Testing FLPQ

Tests are provided as scripts but **do not run automatically**. They require
a compiled FLPQ executable, an MPI runtime, and suitable precomputed OpenMX
fixtures (OpenMX_LPQ.bin and the corresponding OpenMX_LPQ_DM/ directory).
The release archive intentionally does not ship large simulation outputs or
site-specific scheduler logs.

## Methane regression

~~~bash
python3 tests/regression_methane.py --help
python3 tests/regression_methane.py \
  --fixture /path/to/methane-fixture \
  --output-dir /path/to/new-test-results \
  --flpq /path/to/FLPQ/source/flpq
~~~

The output directory must be new; the script copies the fixture data and
checks 1/2/4 MPI rank consistency and low-memory output.

## Stress decomposition

Generate a calculation with Output.Stress.Decomposition on, then:

~~~bash
python3 tests/check_stress_decomposition.py --lpq /path/to/results/LPQ
~~~

This checks stress = stress_S + stress_A, symmetry/antisymmetry and
related scalar quantities.

## HDF5 regression

HDF5=1 and optional NumPy/h5py dependencies are required. The provided
HDF5 scripts were written against a historical fixture layout; they may
require fixture paths matching those referenced in the scripts.

~~~bash
python3 tests/regression_hdf5_direct.py --help
python3 tests/regression_hdf5_f4.py --help
~~~

These scripts exercise direct HDF5 output, quantity selection, precision,
compression, collisions and low-memory consistency, where fixtures are
available. Benchmark helpers are retained separately from correctness
checks; benchmark output must not be interpreted as proof of physical
correctness.

## Quantity-selector unit tests

The selector can be tested without MPI, density-matrix files, or HDF5.
The test compiles the production selector C source with GCC in a
temporary directory under `tests/` and checks all presets, collinear
and noncollinear spin modes, combination rules, and invalid requests:

~~~bash
python3 tests/test_output_quantities.py
~~~

## Baseline release evidence

The release-preparation audit references previously completed P1–P6 tests.
It does not bundle the P1–P6 experimental MPI implementations or their
intermediate data. The HDF5 baseline F3–F5 regressions were performed on
methane and Cr2, including spin/SOC test settings.

For a new release, run at least one clean default HDF5=0 build, one clean
HDF5=1 build, a methane regression, a stress-decomposition case, and
direct-HDF5 tests. Document numerical tolerances and versions used.
