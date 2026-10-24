#!/usr/bin/env python3
"""Validate OpenMX SCF convergence using both energy change and density residual."""

import argparse
import json
from pathlib import Path
import re


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--density-threshold", type=float, default=1e-5)
    parser.add_argument("--energy-threshold", type=float, default=1e-7)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()

    content = args.log.read_text(errors="replace")
    norms = [float(x) for x in re.findall(r"NormRD\s*=\s*([-+0-9.eE]+)", content)]
    changes = [float(x) for x in re.findall(r"dUele\s*=\s*([-+0-9.eE]+)", content)]
    if not norms or not changes or "The calculation was normally finished" not in content:
        raise AssertionError("OpenMX SCF log is incomplete")
    report = {
        "scf_steps": len(norms),
        "last_density_residual": norms[-1],
        "last_energy_change_hartree": changes[-1],
        "density_threshold": args.density_threshold,
        "energy_threshold": args.energy_threshold,
        "converged": norms[-1] < args.density_threshold and abs(changes[-1]) < args.energy_threshold,
    }
    print(json.dumps(report, indent=2))
    if args.report:
        args.report.write_text(json.dumps(report, indent=2)+"\n")
    if not report["converged"]:
        raise AssertionError("SCF did not satisfy both quantitative thresholds")


if __name__ == "__main__":
    main()
