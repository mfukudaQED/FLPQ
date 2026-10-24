#!/usr/bin/env python3
"""Run reproducible OpenMX-to-FLPQ regression tests on an existing fixture.

Requires Python 3.6+, an MPI launcher, and OpenMX LPQ binary fixture files.
The script never modifies the fixture directory.
"""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys


SCALARS = (
    "electron_density",
    "kinetic_energy_density",
    "energy_density_tau",
    "chemical_potential",
)


def sha256(path):
    """Calculate the SHA-256 digest of a file."""
    h = hashlib.sha256()
    with path.open("rb") as stream:
        while True:
            data = stream.read(1024 * 1024)
            if not data:
                break
            h.update(data)
    return h.hexdigest()


def setup_case(fixture, destination, input_name="flpq.inp"):
    """Populate an isolated test directory with immutable fixture copies."""
    if destination.exists():
        raise RuntimeError("Test directory already exists: {}".format(destination))
    destination.mkdir(parents=True)
    shutil.copy2(str(fixture / "OpenMX_LPQ.bin"), str(destination / "OpenMX_LPQ.bin"))
    shutil.copy2(str(fixture / "flpq.inp"), str(destination / input_name))
    shutil.copytree(str(fixture / "OpenMX_LPQ_DM"), str(destination / "OpenMX_LPQ_DM"))


def launch(executable, ranks, directory, input_name, timeout, env, label):
    """Execute FLPQ under mpiexec and record its combined output."""
    command = ["mpiexec", "-n", str(ranks), str(executable), input_name]
    log = directory / (label + ".log")
    with log.open("w") as output:
        try:
            finished = subprocess.run(
                command, cwd=str(directory), stdout=output,
                stderr=subprocess.STDOUT, env=env, timeout=timeout,
            )
            result = finished.returncode
        except subprocess.TimeoutExpired:
            result = 124
    return result, log


def read_numeric_rows(path):
    """Read all finite-valued tabular rows of an LPQ text file."""
    rows = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        numbers = [float(value) for value in stripped.split()]
        if not all(math.isfinite(value) for value in numbers):
            raise AssertionError("Non-finite value in {}".format(path))
        rows.append(numbers)
    return rows


def check_output(directory, grid_count):
    """Validate scalar output counts, finiteness and Cube consistency."""
    output = directory / "LPQ"
    if not output.is_dir():
        raise AssertionError("LPQ directory missing: {}".format(directory))
    for field in SCALARS:
        data = read_numeric_rows(output / ("LPQ." + field + ".dat"))
        if len(data) != grid_count or any(len(row) != 4 for row in data):
            raise AssertionError("Unexpected data grid for {}".format(field))
        cube = (output / ("LPQ." + field + ".cube")).read_text().splitlines()
        atom_count = abs(int(cube[2].split()[0]))
        values = [float(value) for line in cube[6 + atom_count:] for value in line.split()]
        if len(values) != grid_count or not all(math.isfinite(v) for v in values):
            raise AssertionError("Unexpected Cube grid for {}".format(field))
        for row, value in zip(data, values):
            if abs(row[3] - value) > 0.0005 + 0.0005 * abs(row[3]):
                raise AssertionError("Cube/data disagreement for {}".format(field))
    return {path.name: sha256(path) for path in output.iterdir() if path.is_file()}


def compare_grids(reference, candidate, tolerance):
    """Compare numeric LPQ table values across MPI process counts."""
    for field in SCALARS:
        filename = "LPQ." + field + ".dat"
        a = read_numeric_rows(reference / "LPQ" / filename)
        b = read_numeric_rows(candidate / "LPQ" / filename)
        if len(a) != len(b):
            raise AssertionError("MPI grid size mismatch: {}".format(filename))
        for index, (row_a, row_b) in enumerate(zip(a, b)):
            if len(row_a) != len(row_b):
                raise AssertionError("MPI component mismatch: {}".format(filename))
            for x, y in zip(row_a, row_b):
                if abs(x - y) > tolerance * max(1.0, abs(x), abs(y)):
                    raise AssertionError("MPI mismatch at {} point {}".format(filename, index))


def expect_failure(executable, fixture, destination, change, env, timeout):
    """Confirm invalid input fails collectively with a nonzero exit code."""
    setup_case(fixture, destination)
    if change == "missing_dm":
        next((destination / "OpenMX_LPQ_DM").iterdir()).unlink()
    elif change == "short_dm":
        file_path = next((destination / "OpenMX_LPQ_DM").iterdir())
        file_path.write_bytes(file_path.read_bytes()[:8])
    elif change == "short_metadata":
        file_path = destination / "OpenMX_LPQ.bin"
        file_path.write_bytes(file_path.read_bytes()[:8])
    else:
        raise ValueError(change)
    code, log = launch(executable, 2, destination, "flpq.inp", timeout, env, change)
    output = log.read_text(errors="replace")
    if code == 0 or "normally finished" in output:
        raise AssertionError("Invalid input was accepted: {}".format(change))
    if "[FLPQ rank " not in output:
        raise AssertionError("No diagnostic for {}".format(change))
    return code


def main():
    """Run MPI invariance, rerun stability and negative input tests."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--flpq", type=Path, required=True)
    parser.add_argument("--ranks", type=int, nargs="+", default=[1, 2, 4])
    parser.add_argument("--grid-count", type=int, default=729)
    parser.add_argument("--timeout", type=int, default=90)
    parser.add_argument("--tcp-local", action="store_true")
    parser.add_argument("--skip-low-memory", action="store_true", help="Omit slice-mode consistency test")
    args = parser.parse_args()

    fixture = args.fixture.resolve()
    output = args.output_dir.resolve()
    executable = args.flpq.resolve()
    if output.exists():
        parser.error("Refusing to overwrite existing output directory: {}".format(output))
    if not executable.is_file():
        parser.error("Missing FLPQ binary: {}".format(executable))
    output.mkdir(parents=True)

    env = os.environ.copy()
    env["I_MPI_COLL_EXTERNAL"] = "no"
    env["I_MPI_FABRICS"] = "shm:ofi"
    env["OMP_NUM_THREADS"] = "1"
    if args.tcp_local:
        env["I_MPI_OFI_PROVIDER"] = "tcp"
        env["FI_PROVIDER"] = "tcp"

    report = {"ranks": {}, "failure_tests": {}, "tcp_local": args.tcp_local}
    reference = None
    for ranks in args.ranks:
        directory = output / ("ranks_{}".format(ranks))
        setup_case(fixture, directory, input_name="custom_flpq.inp")
        code, log = launch(executable, ranks, directory, "custom_flpq.inp", args.timeout, env, "flpq")
        if code != 0:
            raise AssertionError("MPI {} failed; see {}".format(ranks, log))
        digests = check_output(directory, args.grid_count)
        report["ranks"][str(ranks)] = {"exit_code": code, "digests": digests}
        if reference is None:
            reference = directory
        else:
            compare_grids(reference, directory, 1e-10)
            first = report["ranks"][str(args.ranks[0])]["digests"]
            if digests != first:
                different = sorted(name for name in first if digests.get(name) != first[name])
                raise AssertionError("MPI output byte differences: {}".format(different))

    if args.skip_low_memory:
        report["low_memory_identical"] = None
    else:
        # Verify slice-based low-memory output against the regular 3D output.
        low_memory = output / "low_memory"
        setup_case(fixture, low_memory)
        input_path = low_memory / "flpq.inp"
        input_text = input_path.read_text()
        if "save_memory_3d off" not in input_text:
            raise AssertionError("Fixture must contain save_memory_3d off")
        input_path.write_text(input_text.replace("save_memory_3d off", "save_memory_3d on"))
        code, log = launch(executable, 2, low_memory, "flpq.inp", args.timeout, env, "low_memory")
        if code != 0:
            raise AssertionError("Low-memory test failed; see {}".format(log))
        normal = reference / "LPQ"
        sliced = low_memory / "LPQ"
        for name in sorted(p.name for p in normal.glob("*.dat")):
            if (normal / name).read_bytes() != (sliced / name).read_bytes():
                raise AssertionError("Low-memory output differs: {}".format(name))
        report["low_memory_identical"] = True

    rerun = reference
    before = {p.name: sha256(p) for p in (rerun / "LPQ").iterdir() if p.is_file()}
    code, log = launch(executable, args.ranks[0], rerun, "custom_flpq.inp", args.timeout, env, "rerun")
    if code:
        raise AssertionError("Rerun failed; see {}".format(log))
    after = {p.name: sha256(p) for p in (rerun / "LPQ").iterdir() if p.is_file()}
    if before != after:
        raise AssertionError("Rerun changed output digests")
    report["rerun_identical"] = True

    for invalid in ("missing_dm", "short_dm", "short_metadata"):
        code = expect_failure(executable, fixture, output / invalid, invalid, env, args.timeout)
        report["failure_tests"][invalid] = code

    report["result"] = "PASS"
    (output / "regression_report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print("PASS: MPI {} ranks, {} identical output files, low-memory={}, rerun, 3 invalid inputs".format(args.ranks, len(report["ranks"][str(args.ranks[0])]["digests"]), "SKIPPED" if args.skip_low_memory else "PASS"))
    print("Report: {}".format(output / "regression_report.json"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
