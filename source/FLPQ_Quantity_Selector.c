#define _POSIX_C_SOURCE 200809L
#include "FLPQ_Quantity_Selector.h"

#include <stdio.h>
#include <string.h>
#include <strings.h>

#define FLPQ_SELECTOR_LIMIT 32
#define FLPQ_SPEC_LIMIT 1024

static int selection_error(char *error, size_t error_size,
                           const char *reason, const char *token)
{
  if (error != NULL && error_size > 0)
    snprintf(error, error_size, "%s%s%s",
             reason, token != NULL ? ": " : "",
             token != NULL ? token : "");
  return -1;
}

static int preset_member(const char *token,
                         const FLPQ_Quantity_Descriptor *descriptor)
{
  const char *name = descriptor->name;
  if (!strcasecmp(token, "all_scalar"))
    return descriptor->dimension == 0;
  if (!strcasecmp(token, "all_vector"))
    return descriptor->dimension == 1;
  if (!strcasecmp(token, "all_tensor"))
    return descriptor->dimension == 2;
  if (!strcasecmp(token, "basic"))
    return !strcasecmp(name, "electron_density") ||
           !strcasecmp(name, "spin_angular_momentum") ||
           !strcasecmp(name, "zeta_potential");
  if (!strcasecmp(token, "rcp"))
    return !strcasecmp(name, "kinetic_energy_density") ||
           !strcasecmp(name, "chemical_potential");
  return 0;
}

static int is_preset(const char *token)
{
  return !strcasecmp(token, "all_scalar") ||
         !strcasecmp(token, "all_vector") ||
         !strcasecmp(token, "all_tensor") ||
         !strcasecmp(token, "basic") ||
         !strcasecmp(token, "rcp");
}

int flpq_select_output_quantities(
    const char *request,
    const FLPQ_Quantity_Descriptor *descriptors,
    int count, int spin_mode, int stress_decomposition,
    int *selected, char *error, size_t error_size)
{
  char buffer[FLPQ_SPEC_LIMIT];
  unsigned char seen_explicit[FLPQ_SELECTOR_LIMIT] = {0};
  char *token, *save = NULL;
  int i, count_selected = 0;
  size_t spec_length;

  if (error != NULL && error_size > 0) error[0] = '\0';
  if (!request || !descriptors || !selected || count <= 0 ||
      count > FLPQ_SELECTOR_LIMIT)
    return selection_error(error, error_size, "invalid selection arguments", NULL);

  for (i = 0; i < count; ++i) selected[i] = 0;
  spec_length = strlen(request);
  if (spec_length == 0 || spec_length >= sizeof(buffer))
    return selection_error(error, error_size, "empty or overlong selection", NULL);

  if (!strcasecmp(request, "all")) {
    for (i = 0; i < count; ++i)
      selected[i] =
          (!descriptors[i].noncollinear_only || spin_mode == 3) &&
          (!descriptors[i].optional_stress || stress_decomposition != 0);
    return 0;
  }

  memcpy(buffer, request, spec_length + 1);
  for (token = strtok_r(buffer, ", \t\r\n", &save);
       token != NULL;
       token = strtok_r(NULL, ", \t\r\n", &save)) {
    int match = -1;
    if (!strcasecmp(token, "all"))
      return selection_error(error, error_size, "all cannot be combined", token);

    if (is_preset(token)) {
      for (i = 0; i < count; ++i) {
        if (preset_member(token, &descriptors[i]) &&
            (!descriptors[i].noncollinear_only || spin_mode == 3))
          selected[i] = 1;
      }
    } else {
      for (i = 0; i < count; ++i) {
        if (!strcasecmp(token, descriptors[i].name)) {
          match = i;
          break;
        }
      }
      if (match < 0)
        return selection_error(error, error_size, "unknown quantity or preset", token);
      if (descriptors[match].noncollinear_only && spin_mode != 3)
        return selection_error(error, error_size, "requires noncollinear spin mode", token);
      if (seen_explicit[match])
        return selection_error(error, error_size, "duplicate explicit quantity", token);
      seen_explicit[match] = 1;
      selected[match] = 1;
    }
  }
  for (i = 0; i < count; ++i) count_selected += selected[i] != 0;
  if (count_selected == 0)
    return selection_error(error, error_size, "no quantity selected", NULL);
  return 0;
}
