/* Optional build fallback; legacy FLPQ does not require HDF5. */
#include "flpq.h"
#include "FLPQ_HDF5.h"
#include <stdio.h>
#include <stdlib.h>

void FLPQ_HDF5_write(const char *name, int rank, int index)
{
    int me;
    (void)name; (void)rank; (void)index;
    MPI_Comm_rank(mpi_comm_level1, &me);
    if (me == 0) fprintf(stderr,
        "[FLPQ-HDF5] Output.Format requires rebuilding with make HDF5=1\n");
    MPI_Abort(mpi_comm_level1, EXIT_FAILURE);
}
