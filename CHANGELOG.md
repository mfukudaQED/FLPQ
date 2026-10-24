# Changelog

## v0.2.0 — 2026-10-24

This release adds optional HDF5 output, output validation and conversion
utilities while preserving the established MPI calculation and default
legacy output. The v0.1.0 baseline is the previously published main branch.
Release tags are prepared locally before GitHub publication.

### Added

- Optional direct HDF5 writer, per-quantity HDF5 output selection, and
  converter to legacy dat, Cube and FieldCube.
- HDF5 float64/float32 selection, gzip configuration, chunk/block settings,
  input validation, and collision handling.
- Regression scripts for OpenMX fixtures, stress decomposition, HDF5 output,
  selected quantities and limited-memory execution.
- Updated public build, usage and test documentation.

### Fixed or improved

- OpenMX export and MPI-related error handling and compatibility.
- Stress tensor symmetric/antisymmetric decomposition and associated
  noncollinear/spin accumulation paths.
- Source and output compatibility checks for HDF5 and legacy files.

### Not included

- Experimental P3 asynchronous MPI_Irecv root-HDF5 double buffering.
- Experimental P4 LPT atom reallocation.
- P5/P6 atom-pair, orbital-pair and grid-point distributed task scheduling.

These optimization prototypes did not demonstrate a reliable overall
performance improvement on the available small-molecule cases, and are
kept separately for research rather than incorporated into the release.

### Limitations

- HDF5 writer selection requires a clean rebuild.
- The supplied build flags target an Intel MPI/MKL toolchain; other
  compiler stacks need configuration and validation.
- Fully reproducible test execution requires user-supplied OpenMX fixtures.
