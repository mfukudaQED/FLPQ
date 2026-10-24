#ifndef FLPQ_QUANTITY_SELECTOR_H
#define FLPQ_QUANTITY_SELECTOR_H

#include <stddef.h>

typedef struct {
  const char *name;
  int dimension;
  int noncollinear_only;
  int optional_stress;
} FLPQ_Quantity_Descriptor;

/* Selects a union of explicit quantities and category/purpose presets.
 * Returns 0 on success; otherwise returns -1 and populates error.
 * "all" preserves historical spin/stress defaults and cannot be combined.
 * Spin-incompatible preset members are skipped; explicit requests fail.
 */
int flpq_select_output_quantities(
    const char *request,
    const FLPQ_Quantity_Descriptor *descriptors,
    int count, int spin_mode, int stress_decomposition,
    int *selected, char *error, size_t error_size);

#endif
