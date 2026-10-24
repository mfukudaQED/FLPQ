#ifndef FLPQ_HDF5_WRITER_H
#define FLPQ_HDF5_WRITER_H

/* Collective across the same FLPQ MPI communicator; each call writes exactly
 * one grid slice of a single self-contained per-quantity HDF5 file. */
void FLPQ_HDF5_write(const char *name, int rank, int index);

#endif
