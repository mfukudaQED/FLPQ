#!/usr/bin/env python3
"""Convert self-contained FLPQ-HDF5 v1 fields to cube, FieldCube v1, or legacy dat.

This is an independent post-processing CLI; it does not load OpenMX/FLPQ inputs.
The HDF5 datasets are read one first-axis grid slab at a time.
"""

import argparse
import json
import os
from pathlib import Path
import sys
import tempfile

import h5py
import numpy as np


BOHR_TO_ANGSTROM = 0.529177210903
AXES = ("x", "y", "z")
TENSOR_COMPONENTS = ("xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz")
REQUIRED_ATTRS = (
    "format_id", "schema_version", "producer", "quantity", "coordinate_unit",
    "value_unit", "field_rank", "precision", "component_basis",
    "component_order", "data_order",
)


class FormatError(ValueError):
    """HDF5 data does not satisfy the FLPQ v1 file contract."""


def string_attr(value):
    """Normalize HDF5 byte/string scalar attributes to Python strings."""
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, str):
        return value
    raise FormatError("Expected a string attribute; got {!r}".format(value))


def require_array(file, name, shape, kind=None):
    if name not in file or not isinstance(file[name], h5py.Dataset):
        raise FormatError("Missing HDF5 dataset: " + name)
    ds = file[name]
    if len(ds.shape) != len(shape) or any(
        target is not None and got != target for got, target in zip(ds.shape, shape)
    ):
        raise FormatError("{} has shape {}; expected {}".format(name, ds.shape, shape))
    if kind and ds.dtype.kind not in kind:
        raise FormatError("{} has unsupported dtype {}".format(name, ds.dtype))
    return ds


class FLPQFile:
    """Validated metadata and lazily accessible HDF5 datasets."""

    def __init__(self, path):
        self.path = Path(path)
        self.file = h5py.File(str(self.path), "r")
        try:
            self._validate()
        except BaseException:
            self.file.close()
            raise

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback):
        self.file.close()

    def _validate(self):
        h = self.file
        for name in REQUIRED_ATTRS:
            if name not in h.attrs:
                raise FormatError("Missing required attribute: " + name)
        if string_attr(h.attrs["format_id"]) != "FLPQ-HDF5":
            raise FormatError("Not an FLPQ-HDF5 file")
        if int(h.attrs["schema_version"]) != 1:
            raise FormatError("Unsupported FLPQ-HDF5 schema version")
        if string_attr(h.attrs["coordinate_unit"]).lower() != "bohr":
            raise FormatError("v1 expects coordinates in bohr")
        if string_attr(h.attrs["component_basis"]) != "cartesian":
            raise FormatError("v1 expects Cartesian field components")
        if string_attr(h.attrs["data_order"]) != "component,n1,n2,n3":
            raise FormatError("Unsupported HDF5 data order")
        self.quantity = string_attr(h.attrs["quantity"])
        self.value_unit = string_attr(h.attrs["value_unit"])
        self.field_rank = int(h.attrs["field_rank"])
        self.precision = string_attr(h.attrs["precision"])
        self.component_order = string_attr(h.attrs["component_order"])
        if self.precision not in ("float32", "float64"):
            raise FormatError("Unsupported precision {}".format(self.precision))
        if self.quantity == "stress_diag":
            if self.field_rank != 2:
                raise FormatError("stress_diag must have rank 2")
        elif self.field_rank not in (0, 1, 2):
            raise FormatError("Invalid field rank")

        shape = np.asarray(require_array(h, "/grid/shape", (3,), "iu")[:], dtype=np.int64)
        if np.any(shape <= 0):
            raise FormatError("Grid dimensions must be positive")
        self.shape = tuple(map(int, shape))
        self.origin = np.asarray(require_array(h, "/grid/origin", (3,), "f")[:], dtype=np.float64)
        self.steps = np.asarray(require_array(h, "/grid/step_vectors", (3, 3), "f")[:], dtype=np.float64)
        self.lattice = np.asarray(require_array(h, "/grid/lattice_vectors", (3, 3), "f")[:], dtype=np.float64)
        self.periodic = np.asarray(require_array(h, "/grid/periodic", (3,), "iub")[:], dtype=int)
        if not all(np.all(np.isfinite(a)) for a in (self.origin, self.steps, self.lattice)):
            raise FormatError("Non-finite grid geometry")
        if np.any((self.periodic != 0) & (self.periodic != 1)):
            raise FormatError("periodic entries must be 0 or 1")
        if np.linalg.matrix_rank(self.steps) != 3:
            raise FormatError("Grid step vectors must span 3D space")

        atoms = require_array(h, "/atoms/atomic_numbers", (None,), "iu")
        count = atoms.shape[0]
        require_array(h, "/atoms/positions", (count, 3), "f")
        require_array(h, "/atoms/cube_charges", (count,), "f")
        self.natoms = count
        self.atomic_numbers = np.asarray(atoms[:], dtype=np.int32)
        self.positions = np.asarray(h["/atoms/positions"][:], dtype=np.float64)
        self.cube_charges = np.asarray(h["/atoms/cube_charges"][:], dtype=np.float64)
        if (np.any(self.atomic_numbers < 0) or np.any(self.atomic_numbers > 118)
                or not np.all(np.isfinite(self.positions))
                or not np.all(np.isfinite(self.cube_charges))):
            raise FormatError("Invalid atomic metadata")

        dtype_itemsize = {"float32": 4, "float64": 8}[self.precision]
        if self.quantity == "stress_diag":
            self.values = None
            self.eigenvalues = require_array(h, "/field/eigenvalues", (3,) + self.shape, "f")
            self.eigenvectors = require_array(h, "/field/eigenvectors", (3, 3) + self.shape, "f")
            datasets = [self.eigenvalues, self.eigenvectors]
            self.ncomp = 12
        else:
            self.ncomp = {0: 1, 1: 3, 2: 9}[self.field_rank]
            expected = {0: "scalar", 1: "x y z", 2: " ".join(TENSOR_COMPONENTS)}[self.field_rank]
            if self.component_order != expected:
                raise FormatError("Component order mismatch: expected {!r}".format(expected))
            self.values = require_array(h, "/field/values", (self.ncomp,) + self.shape, "f")
            datasets = [self.values]
            self.eigenvalues = self.eigenvectors = None
        for ds in datasets:
            if ds.dtype.itemsize != dtype_itemsize:
                raise FormatError("{} dtype differs from precision {}".format(ds.name, self.precision))

    def scalar_index(self, component, eigenvalue):
        """Resolve exactly one scalar component for cube generation."""
        if self.quantity == "stress_diag":
            if component is not None:
                raise FormatError("Use --eigenvalue 1/2/3 for stress_diag")
            if eigenvalue not in (1, 2, 3):
                raise FormatError("stress_diag requires --eigenvalue 1, 2 or 3")
            return ("eigenvalues", eigenvalue - 1)
        if eigenvalue is not None:
            raise FormatError("--eigenvalue is only valid for stress_diag")
        if self.field_rank == 0:
            if component not in (None, "scalar"):
                raise FormatError("Scalar field has no component {!r}".format(component))
            return ("values", 0)
        choices = AXES if self.field_rank == 1 else TENSOR_COMPONENTS
        if component not in choices:
            raise FormatError("Select --component from: " + ", ".join(choices))
        return ("values", choices.index(component))

    def slab(self, index):
        """Return [component, j, k] for a constant first-grid-axis index."""
        if self.quantity == "stress_diag":
            ev = np.asarray(self.eigenvalues[:, index, :, :], dtype=np.float64)
            vec = np.asarray(self.eigenvectors[:, :, index, :, :], dtype=np.float64)
            slab = np.concatenate((ev, vec.reshape(9, *self.shape[1:])), axis=0)
        else:
            slab = np.asarray(self.values[:, index, :, :], dtype=np.float64)
        if not np.all(np.isfinite(slab)):
            raise FormatError("Non-finite field values at grid plane {}".format(index))
        return slab

    def info(self):
        return {
            "format_id": "FLPQ-HDF5", "schema_version": 1,
            "quantity": self.quantity, "field_rank": self.field_rank,
            "components": self.ncomp, "shape": list(self.shape),
            "precision": self.precision, "coordinate_unit": "bohr",
            "value_unit": self.value_unit, "atoms": self.natoms,
            "component_order": self.component_order,
            "datasets": (["/field/eigenvalues", "/field/eigenvectors"]
                         if self.quantity == "stress_diag" else ["/field/values"]),
        }


def atomic_output(destination, producer, overwrite=False):
    """Write a completed output using a temporary file in its own directory."""
    output = Path(destination)
    if output.is_dir():
        raise FormatError("Output path is a directory")
    if not output.parent.is_dir():
        raise FormatError("Output directory does not exist: {}".format(output.parent))
    if output.exists() and not overwrite:
        raise FileExistsError("Output already exists (use --overwrite): {}".format(output))
    name = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", newline="\n",
            prefix=".flpq-convert-", suffix=".partial",
            dir=str(output.parent), delete=False
        ) as fp:
            name = fp.name
            producer(fp)
            fp.flush()
            os.fsync(fp.fileno())
        if overwrite:
            os.replace(name, str(output))
        else:
            # Atomic no-clobber publish on the same user-controlled filesystem.
            os.link(name, str(output))
            os.unlink(name)
        name = None
    finally:
        if name and os.path.exists(name):
            os.unlink(name)


def write_cube_field(fp, data, selector, digits):
    """Emit scalar, 6 numbers per line, k fastest (j rows separated)."""
    where, icomp = selector
    n1, n2, n3 = data.shape
    fmt = "{:." + str(digits) + "E}"
    for i in range(n1):
        slab = data.slab(i)
        for j in range(n2):
            row = slab[icomp, j, :]
            for k in range(0, n3, 6):
                fp.write(" ".join(fmt.format(float(v)) for v in row[k:k+6]) + "\n")


def export_cube(data, fp, selector, digits):
    fp.write("FLPQ-HDF5 v1 cube export: {}\n".format(data.quantity))
    fp.write("Converted from HDF5, coordinates in bohr, field in {}\n".format(data.value_unit))
    fp.write("{:5d} {: .10f} {: .10f} {: .10f}\n".format(data.natoms, *data.origin))
    for n, vec in zip(data.shape, data.steps):
        fp.write("{:5d} {: .10f} {: .10f} {: .10f}\n".format(n, *vec))
    for z, charge, r in zip(data.atomic_numbers, data.cube_charges, data.positions):
        fp.write("{:5d} {: .10f} {: .10f} {: .10f} {: .10f}\n".format(int(z), charge, *r))
    write_cube_field(fp, data, selector, digits)


def export_fieldcube(data, fp, digits):
    if data.quantity == "stress_diag":
        raise FormatError("stress_diag is not an ordinary tensor FieldCube; use cube --eigenvalue or dat")
    meta = {
        "name": data.quantity,
        "fieldType": {0: "scalar", 1: "vector", 2: "tensor"}[data.field_rank],
        "rank": data.field_rank,
        "shape": {0: [], 1: [3], 2: [3, 3]}[data.field_rank],
        "components": data.ncomp,
        "componentOrder": data.component_order,
        "basis": "cartesian",
        "dataOrder": "i-j-k-kfastest",
        "coordinateUnit": "Ang",
        "valueUnit": data.value_unit,
        "periodic": bool(np.any(data.periodic)),
        "symmetry": "none",
    }
    fp.write("PQViewer FieldCube v1\n")
    fp.write(json.dumps(meta, ensure_ascii=False, separators=(",", ":")) + "\n")
    origin = data.origin * BOHR_TO_ANGSTROM
    fp.write("{:d} {:.14E} {:.14E} {:.14E}\n".format(data.natoms, *origin))
    for n, vec in zip(data.shape, data.steps * BOHR_TO_ANGSTROM):
        fp.write("{:d} {:.14E} {:.14E} {:.14E}\n".format(n, *vec))
    for z, charge, r in zip(data.atomic_numbers, data.cube_charges,
                             data.positions * BOHR_TO_ANGSTROM):
        fp.write("{} {:.14E} {:.14E} {:.14E} {:.14E}\n".format(int(z), charge, *r))
    fmt = "{:." + str(digits) + "E}"
    for i in range(data.shape[0]):
        slab = data.slab(i)
        for point in np.moveaxis(slab, 0, -1).reshape(-1, data.ncomp):
            fp.write(" ".join(fmt.format(float(x)) for x in point) + "\n")


def export_dat(data, fp, digits, legacy_format):
    n1, n2, n3 = data.shape
    fp.write("# Ngrid {} {} {}\n".format(n1, n2, n3))
    if data.quantity == "stress_diag":
        note = "eigenvalues(3), eigenvectors(9)"
    else:
        note = {0: "scalar", 1: "vector (x,y,z)",
                2: "tensor (xx,yx,zx,xy,yy,zy,xz,yz,zz)"}[data.field_rank]
    fp.write("# x, y, z, {} (unit : [bohr], [{}])\n".format(note, data.value_unit))
    fval = ("{:18.8E}" if legacy_format else "{:." + str(digits) + "E}")
    for i in range(n1):
        slab = data.slab(i)
        for j in range(n2):
            for k in range(n3):
                r = data.origin + i*data.steps[0] + j*data.steps[1] + k*data.steps[2]
                if legacy_format:
                    coords = "".join("{:14.4E}".format(float(x)) for x in r)
                else:
                    coords = " ".join("{:.14E}".format(float(x)) for x in r) + " "
                point = slab[:, j, k]
                if data.field_rank == 2:
                    if data.quantity == "stress_diag":
                        # Eigenvalue(3), then eigenvector matrix column-major.
                        matrix = point[3:].reshape(3, 3)
                        point = np.concatenate([point[:3], matrix.T.reshape(9)])
                    else:
                        # HDF5 tensor is row-major; legacy FLPQ dat is column-major.
                        point = point.reshape(3, 3).T.reshape(9)
                values_text = ("".join(fval.format(float(v)) for v in point)
                               if legacy_format else " ".join(fval.format(float(v)) for v in point))
                fp.write(coords + values_text + "\n")
            fp.write("\n")


def build_parser():
    parser = argparse.ArgumentParser(
        description="Inspect/convert independent FLPQ-HDF5 v1 field files"
    )
    sub = parser.add_subparsers(dest="command", required=True)
    view = sub.add_parser("inspect", help="Validate and print HDF5 schema information")
    view.add_argument("input", type=Path)
    for name in ("cube", "fieldcube", "dat"):
        p = sub.add_parser(name, help="Convert FLPQ-HDF5 to {}".format(name))
        p.add_argument("input", type=Path)
        p.add_argument("--output", type=Path, required=True)
        p.add_argument("--overwrite", action="store_true")
        p.add_argument("--digits", type=int, default=14, help="Decimal digits after E mantissa, 3..17")
        if name == "cube":
            p.add_argument("--component", help="Vector x/y/z or tensor xx/xy/.../zz")
            p.add_argument("--eigenvalue", type=int, choices=(1, 2, 3))
        if name == "dat":
            p.add_argument("--legacy-format", action="store_true",
                           help="Use legacy FLPQ fixed-width numeric output")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command != "inspect" and not (3 <= args.digits <= 17):
            raise FormatError("--digits must be in [3, 17]")
        with FLPQFile(args.input) as data:
            if args.command == "inspect":
                print(json.dumps(data.info(), sort_keys=True, indent=2))
                return 0
            if args.output.resolve() == args.input.resolve():
                raise FormatError("Refusing to overwrite the input HDF5 file")
            if args.command == "cube":
                selector = data.scalar_index(args.component, args.eigenvalue)
                output = lambda fp: export_cube(data, fp, selector, args.digits)
            elif args.command == "fieldcube":
                if data.quantity == "stress_diag":
                    raise FormatError("stress_diag FieldCube conversion unsupported; use cube --eigenvalue or dat")
                output = lambda fp: export_fieldcube(data, fp, args.digits)
            else:
                output = lambda fp: export_dat(data, fp, args.digits, args.legacy_format)
            atomic_output(args.output, output, args.overwrite)
            print("Wrote {}".format(args.output))
            return 0
    except (OSError, ValueError, TypeError, KeyError) as exc:
        print("FLPQ conversion error: {}".format(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
