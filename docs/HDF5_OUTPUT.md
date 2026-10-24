# Optional HDF5 output

Direct HDF5 writing is available when FLPQ is compiled using HDF5=1.
HDF5=0 and unspecified Output.Format preserve the legacy behavior.

## Example

~~~text
Output.Format hdf5
Output.Quantities electron_density,charge_current,stress
Output.HDF5.Precision float64
Output.HDF5.Compression gzip
Output.HDF5.GzipLevel 4
Output.HDF5.OnExisting error
~~~

HDF5 is self-contained, one file per selected physical quantity. Each
file stores scalar/vector/tensor field datasets, grid origin and steps,
atom positions, labels, field unit metadata and component ordering.

## Settings

| Input name | Values | Default |
| --- | --- | --- |
| Output.Format | legacy, hdf5, both | legacy |
| Output.Quantities | all, all_scalar, all_vector, all_tensor, basic, rcp, or comma/space-separated quantities/presets | all |
| Output.HDF5.Precision | float64, float32 | float64 |
| Output.HDF5.Compression | gzip, none | gzip |
| Output.HDF5.GzipLevel | 0–9 | 4 |
| Output.HDF5.OnExisting | error, replace | error |
| Output.HDF5.ChunkY | 1–128 (bounded by grid size) | 16 |
| Output.HDF5.ChunkZ | 1–512 (bounded by grid size) | 128 |
| Output.HDF5.BlockPoints | 16–8192 | 512 |
| Output.HDF5.LogPerformance | on, off | off |

See the **[complete quantity and preset reference](OUTPUT_QUANTITIES.md)**
for all 14 names, scalar/vector/tensor groupings, the `basic` and `rcp`
presets, combinations, and spin-mode restrictions. The default `all`
retains historical output behavior. Presets skip spin-incompatible members,
whereas explicitly requesting an incompatible individual quantity is an
error. Duplicate explicit quantity names and `all` mixed with other tokens
remain errors.

The stress_diag data contains three eigenvalues and nine eigenvector
components. It is not simply a nine-component stress tensor. Eigenvector
signs may vary between MPI layouts due to roundoff and near degeneracy.

## Precision and compression

Float64 is the preferred choice for numerical validation; float32 rounds
the calculated values and can reduce storage. Gzip is lossless; higher
levels may consume more CPU time. The HDF5 writer is root-based; parallel
HDF5 support is not required.

## Existing files and partial writes

By default, Output.HDF5.OnExisting error refuses to overwrite finished
HDF5 files. The replace option writes partial files before replacing
the corresponding completed output. Existing partial files are not
silently deleted. File replacement is **per quantity**, not an atomic
transaction over multiple quantities, so an interrupted run can leave
different output files from different attempts.

OnExisting does not affect legacy dat/Cube collision behavior. Do not
run concurrent FLPQ calculations into the same output directory.

## Conversion

Install the optional Python dependencies:

~~~bash
python3 -m pip install -r tools/requirements-hdf5.txt
python3 tools/flpq_convert.py --help
python3 tools/flpq_convert.py inspect path/to/LPQ.electron_density.h5
~~~

The converter reads HDF5 in bounded slabs and can export Cube, FieldCube,
or legacy dat. For vectors and tensors, the requested component and
ordering must be checked against the metadata. Preserve the original
HDF5 file if you need full precision.
