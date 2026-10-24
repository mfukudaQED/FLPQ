# Selectable physical quantities

Use the input keyword `Output.Quantities` in `flpq.inp` to control
which physical quantities FLPQ evaluates and writes. The selection applies
to **both legacy dat/Cube output and HDF5 output**, independent of
`Output.Format`. The default is `Output.Quantities all`.

## Supported quantities

| Category | Quantity name | Notes |
| --- | --- | --- |
| Scalar | `electron_density` | Electron density |
| Scalar | `kinetic_energy_density` | Kinetic energy density |
| Scalar | `energy_density_tau` | Energy density from stress trace |
| Scalar | `chemical_potential` | Local chemical potential |
| Scalar | `zeta_potential` | Noncollinear spin only (`SpinP_switch == 3`) |
| Vector | `charge_current` | Charge-current density |
| Vector | `kinetic_momentum` | Kinetic-momentum density |
| Vector | `spin_vorticity` | Spin vorticity |
| Vector | `spin_angular_momentum` | Spin angular momentum density |
| Vector | `zeta_force` | Noncollinear spin only (`SpinP_switch == 3`) |
| Tensor | `stress_diag` | Three eigenvalues and associated eigenvectors |
| Tensor | `stress` | Full electronic stress tensor |
| Tensor | `stress_S` | Symmetric part of the stress tensor |
| Tensor | `stress_A` | Antisymmetric part of the stress tensor |

These are the 14 registered outputs, not all internal numerical evaluation
functions implemented by FLPQ. `stress_diag` is classified as a tensor
output although its representation comprises eigenvalues and eigenvectors,
rather than a conventional 3-by-3 stress matrix.

## Selection presets

| Preset | Quantities selected |
| --- | --- |
| `all` (default) | Historical output set for the active spin mode; includes `stress`, `stress_S`, and `stress_A` only when `Output.Stress.Decomposition on` |
| `all_scalar` | All scalar quantities listed above |
| `all_vector` | All vector quantities listed above |
| `all_tensor` | All four tensor outputs listed above |
| `basic` | `electron_density`, `spin_angular_momentum`, `zeta_potential` |
| `rcp` | `kinetic_energy_density`, `chemical_potential` |

The presets `all_scalar`, `all_vector`, and `basic` omit
noncollinear-only quantities when the input does not use noncollinear spin:
for example, `basic` then produces only `electron_density` and
`spin_angular_momentum`. In noncollinear calculations, `basic` produces
all three specified quantities. `all_tensor` includes `stress`,
`stress_S`, and `stress_A` even if
`Output.Stress.Decomposition` is `off`.

Individual names and presets may be combined with commas or whitespace.
Their union is selected, with overlaps between presets deduplicated.
Explicitly repeating the same individual name remains an error.
`all` is a **standalone** selection and cannot be mixed with other names.
Requesting a noncollinear-only quantity *explicitly* in an incompatible
spin mode is an error, as in the original quantity selector.

### Examples

```text
Output.Quantities basic
Output.Quantities rcp
Output.Quantities all_scalar
Output.Quantities all_vector,all_tensor
Output.Quantities electron_density,chemical_potential
Output.Quantities basic,rcp
```

Use **only one** `Output.Quantities` line per input file; the examples
above show alternatives, not a single valid input.

The selection changes which derived local physical quantities are
evaluated and written. It does not yet disable preparation of shared
orbital/density/derivative grids in the upstream FLPQ calculation.

For output encoding, compression and file-collision behavior, see
[HDF5 output](HDF5_OUTPUT.md).
