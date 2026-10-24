#!/usr/bin/env python3
"""Verify stress = symmetric + antisymmetric for every FLPQ grid point."""

import argparse
import math
from pathlib import Path


def rows(path, components):
    """Parse FLPQ's three-position-column text format."""
    data = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        values = [float(value) for value in line.split()]
        if len(values) != components + 3 or not all(math.isfinite(v) for v in values):
            raise AssertionError("Invalid row in {}: {}".format(path.name, line[:100]))
        data.append(values)
    if not data:
        raise AssertionError("Empty output: {}".format(path))
    return data


def tensor(values):
    """Output ordering is xx,yx,zx,xy,yy,zy,xz,yz,zz."""
    flat = values[3:]
    return [[flat[3 * j + i] for j in range(3)] for i in range(3)]


def near(a, b, atol=4e-8, rtol=3e-7):
    return abs(a - b) <= atol + rtol * max(abs(a), abs(b))


def check(folder, require_nonzero):
    full = rows(folder / "LPQ.stress.dat", 9)
    symmetric = rows(folder / "LPQ.stress_S.dat", 9)
    antisymmetric = rows(folder / "LPQ.stress_A.dat", 9)
    energies = rows(folder / "LPQ.energy_density_tau.dat", 1)
    potentials = rows(folder / "LPQ.chemical_potential.dat", 1)
    densities = rows(folder / "LPQ.electron_density.dat", 1)
    if len({len(x) for x in (full,symmetric,antisymmetric,energies,potentials,densities)})!=1:
        raise AssertionError("Grid lengths disagree")
    max_antisym = 0.0
    max_residual = 0.0
    for k, (fr,sr,ar,er,mr,dr) in enumerate(zip(full,symmetric,antisymmetric,energies,potentials,densities)):
        for other in (sr,ar,er,mr,dr):
            for x,y in zip(fr[:3], other[:3]):
                if not near(x,y,1e-12,0):
                    raise AssertionError("Grid coordinates differ at {}".format(k))
        f,s,a = tensor(fr),tensor(sr),tensor(ar)
        for i in range(3):
            for j in range(3):
                max_antisym=max(max_antisym,abs(a[i][j]))
                max_residual=max(max_residual,abs(f[i][j]-s[i][j]-a[i][j]))
                if not near(f[i][j],s[i][j]+a[i][j]):
                    raise AssertionError("Sigma != S + A, row {} component {} {}".format(k,i,j))
                if not near(s[i][j],s[j][i]):
                    raise AssertionError("Symmetry violation, row {}".format(k))
                if not near(a[i][j],-a[j][i]):
                    raise AssertionError("Antisymmetry violation, row {}".format(k))
            if not near(a[i][i],0.0,1e-12,0):
                raise AssertionError("Nonzero antisymmetric diagonal at {}".format(k))
        expected_energy=0.5*sum(f[i][i] for i in range(3))
        if not near(er[3],expected_energy):
            raise AssertionError("Energy != 0.5 trace at {}".format(k))
        density=dr[3]
        expected_mu=0.0 if abs(density)<1e-16 else expected_energy/density
        if abs(density)>=1e-12 and not near(mr[3],expected_mu,4e-7,6e-7):
            raise AssertionError("RCP != energy / density at {}".format(k))
    if require_nonzero and max_antisym<=1e-9:
        raise AssertionError("Antisymmetric stress has no measurable nonzero component")
    print("PASS: {} points, max antisymmetric={:.7g}, max sum residual={:.7g}".format(
        len(full),max_antisym,max_residual))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lpq",type=Path,required=True)
    parser.add_argument("--require-nonzero-antisymmetric",action="store_true")
    args=parser.parse_args()
    check(args.lpq,args.require_nonzero_antisymmetric)


if __name__=="__main__":
    main()
