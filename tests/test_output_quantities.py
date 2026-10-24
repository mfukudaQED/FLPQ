#!/usr/bin/env python3
"""Compile and exercise the C quantity selector without MPI or OpenMX fixtures."""

import ctypes
from pathlib import Path
import subprocess
import tempfile


class Descriptor(ctypes.Structure):
    """Match FLPQ_Quantity_Descriptor from the C header."""

    _fields_ = [
        ("name", ctypes.c_char_p),
        ("dimension", ctypes.c_int),
        ("noncollinear_only", ctypes.c_int),
        ("optional_stress", ctypes.c_int),
    ]


SCALARS = [
    "electron_density",
    "zeta_potential",
    "kinetic_energy_density",
    "energy_density_tau",
    "chemical_potential",
]
VECTORS = [
    "charge_current",
    "kinetic_momentum",
    "spin_vorticity",
    "spin_angular_momentum",
    "zeta_force",
]
TENSORS = ["stress_diag", "stress", "stress_S", "stress_A"]
NAMES = SCALARS + VECTORS + TENSORS
DESCRIPTORS = (Descriptor * len(NAMES))(
    *[
        Descriptor(
            name.encode(),
            0 if name in SCALARS else (1 if name in VECTORS else 2),
            int(name in {"zeta_potential", "zeta_force"}),
            int(name in {"stress", "stress_S", "stress_A"}),
        )
        for name in NAMES
    ]
)


def load_selector(repository: Path, build_dir: Path):
    """Build the production C selector as a standalone shared library."""
    source = repository / "source" / "FLPQ_Quantity_Selector.c"
    binary = build_dir / "libflpq_selector.so"
    subprocess.run(
        [
            "gcc", "-std=c11", "-Wall", "-Wextra", "-Werror", "-O2",
            "-fPIC", "-shared", str(source), "-o", str(binary),
        ],
        check=True,
    )
    method = ctypes.CDLL(str(binary)).flpq_select_output_quantities
    method.argtypes = [
        ctypes.c_char_p,
        ctypes.POINTER(Descriptor),
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_char_p,
        ctypes.c_size_t,
    ]
    method.restype = ctypes.c_int
    return method


def select(method, spec: str, spin: int, decomposition: int = 0):
    """Return the selected names or the C parser's reported error."""
    chosen = (ctypes.c_int * len(NAMES))()
    error = ctypes.create_string_buffer(256)
    status = method(
        spec.encode(),
        DESCRIPTORS,
        len(NAMES),
        spin,
        decomposition,
        chosen,
        error,
        len(error),
    )
    if status != 0:
        return None, error.value.decode()
    return {name for i, name in enumerate(NAMES) if chosen[i]}, None


def verify(method, spec, spin, expected, decomposition=0):
    """Assert exact quantity selection for one input and spin mode."""
    actual, error = select(method, spec, spin, decomposition)
    assert error is None, (spec, spin, error)
    assert actual == set(expected), (spec, spin, actual, expected)


def run_tests(method):
    """Check presets, combinations, compatibility, and error handling."""
    baseline = set(SCALARS + VECTORS + ["stress_diag"])
    noncollinear = set(baseline)
    baseline -= {"zeta_potential", "zeta_force"}

    verify(method, "all", 0, baseline)
    verify(method, "all", 3, noncollinear)
    verify(method, "all", 0, baseline | set(TENSORS), decomposition=1)
    verify(method, "all_scalar", 0, set(SCALARS) - {"zeta_potential"})
    verify(method, "all_scalar", 3, SCALARS)
    verify(method, "all_vector", 0, set(VECTORS) - {"zeta_force"})
    verify(method, "all_vector", 3, VECTORS)
    verify(method, "all_tensor", 0, TENSORS)
    verify(method, "all_tensor", 3, TENSORS)
    verify(method, "basic", 0, ["electron_density", "spin_angular_momentum"])
    verify(
        method, "basic", 3,
        ["electron_density", "spin_angular_momentum", "zeta_potential"],
    )
    verify(method, "rcp", 0, ["kinetic_energy_density", "chemical_potential"])
    verify(method, "rcp", 3, ["kinetic_energy_density", "chemical_potential"])
    verify(method, "all_scalar,all_vector", 0, baseline - {"stress_diag"})
    verify(method, "basic, rcp", 3, [
        "electron_density", "spin_angular_momentum", "zeta_potential",
        "kinetic_energy_density", "chemical_potential",
    ])
    verify(method, "all_tensor,chemical_potential", 0, TENSORS + ["chemical_potential"])
    verify(method, "RCP", 0, ["kinetic_energy_density", "chemical_potential"])
    verify(method, "all_scalar,rcp", 0, set(SCALARS) - {"zeta_potential"})
    verify(method, "electron_density,charge_current", 0, [
        "electron_density", "charge_current",
    ])
    for value in [
        "", "  ", "all,rcp", "electron_density,electron_density",
        "not_a_quantity", "zeta_potential", "zeta_force",
    ]:
        selected, error = select(method, value, 0)
        assert selected is None and error, (value, selected, error)
    print("PASS: C selector presets, spin modes, union and invalid inputs")


def main():
    """Run tests in an isolated build directory under the project."""
    repo = Path(__file__).resolve().parents[1]
    with tempfile.TemporaryDirectory(prefix=".flpq-selector-", dir=repo / "tests") as directory:
        method = load_selector(repo, Path(directory))
        run_tests(method)


if __name__ == "__main__":
    main()
