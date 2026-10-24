#!/usr/bin/env python3
"""F4 real-FLPQ acceptance: selected quantities, encodings and collision safety."""
import argparse
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import h5py
import numpy as np

from regression_hdf5_direct import sha256, rows, field_array
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"tools"))
from flpq_convert import FLPQFile

DEFAULT_NAMES_METHANE={
    "electron_density","kinetic_energy_density","energy_density_tau",
    "chemical_potential","charge_current","kinetic_momentum","spin_vorticity",
    "spin_angular_momentum","stress_diag"
}
DEFAULT_NAMES_CR2=DEFAULT_NAMES_METHANE|{"zeta_potential","zeta_force"}
STRESS_NAMES={"stress","stress_S","stress_A"}


def prepare_case(fixture,target,append="",low=False,original=True):
    if target.exists():raise AssertionError("Existing test case: "+str(target))
    target.mkdir(parents=True)
    shutil.copy2(fixture/"OpenMX_LPQ.bin",target/"OpenMX_LPQ.bin")
    shutil.copytree(fixture/"OpenMX_LPQ_DM",target/"OpenMX_LPQ_DM")
    s=(fixture/"flpq.inp").read_text()
    assert "Output.Format" not in s and "Output.Quantities" not in s
    if low:
        assert "save_memory_3d off" in s
        s=s.replace("save_memory_3d off","save_memory_3d on")
    (target/"flpq.inp").write_text(s+("\n"+append+"\n" if append else ""))


def append_case(case,body):
    input=case/"flpq.inp"
    input.write_text(input.read_text()+"\n"+body+"\n")


def run_case(exe,case,ranks=2,expect_ok=True):
    command=["mpiexec","-n",str(ranks),str(exe),"flpq.inp"]
    with (case/"run.log").open("w") as log:
        p=subprocess.run(command,cwd=str(case),stdout=log,
                         stderr=subprocess.STDOUT,timeout=65)
    text=(case/"run.log").read_text(errors="replace")
    if expect_ok and p.returncode!=0:
        raise AssertionError("unexpected failure {}: {}".format(case,text[-3000:]))
    if not expect_ok and p.returncode==0:
        raise AssertionError("invalid input/collision unexpectedly passed "+str(case))
    return text


def outputs(case):
    p=case/"LPQ"
    return {f.name for f in p.iterdir() if f.is_file()} if p.exists() else set()


def check_subset(case,expected,reference,precision="float64",compression="gzip",
                 gzip_level=4,require_legacy=False):
    found={name[len("LPQ."):-3] for name in outputs(case)
           if name.startswith("LPQ.") and name.endswith(".h5")}
    assert found==set(expected),(case,"HDF5",found,expected)
    for name in expected:
        h5=case/"LPQ"/("LPQ."+name+".h5")
        with FLPQFile(h5) as field:
            assert field.quantity==name
            assert field.precision==precision
            dset=field.file["/field/eigenvalues" if name=="stress_diag" else "/field/values"]
            assert dset.compression==(None if compression=="none" else "gzip")
            if compression=="gzip":assert dset.compression_opts==gzip_level
            stored_compression=field.file.attrs["compression"]
            if isinstance(stored_compression,bytes):
                stored_compression=stored_compression.decode("utf-8")
            assert stored_compression==compression,(stored_compression,compression)
            assert field.file.attrs["gzip_level"]==(gzip_level if compression=="gzip" else 0)
            p=reference/("LPQ."+name+".dat")
            original=rows(p)
            result=field_array(field)
            assert result.shape==original[:,3:].shape,(name,result.shape,original.shape)
            if precision=="float32":
                np.testing.assert_allclose(result,original[:,3:],rtol=3e-6,atol=3e-7,
                                           err_msg=name+" float32")
            else:
                np.testing.assert_allclose(result,original[:,3:],rtol=6e-9,atol=1e-12,
                                           err_msg=name+" float64")
        if require_legacy:
            assert sha256(case/"LPQ"/("LPQ."+name+".dat"))==sha256(reference/("LPQ."+name+".dat")),name
    if not require_legacy:
        assert not any(f.endswith(".dat") or f.endswith(".cube") for f in outputs(case)),case
    else:
        legacy_dat={x for x in outputs(case) if x.endswith(".dat")}
        assert legacy_dat=={"LPQ."+name+".dat" for name in expected},legacy_dat
    assert not list((case/"LPQ").glob("*.partial.h5"))


def compare_selected(a,b,names):
    for name in names:
        with FLPQFile(a/"LPQ"/("LPQ."+name+".h5")) as x, FLPQFile(b/"LPQ"/("LPQ."+name+".h5")) as y:
            assert x.shape==y.shape
            np.testing.assert_allclose(field_array(x),field_array(y),rtol=3e-13,atol=1e-13)
            np.testing.assert_allclose(x.origin,y.origin,rtol=0,atol=1e-13)
            np.testing.assert_allclose(x.steps,y.steps,rtol=0,atol=1e-13)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--base",required=True,type=Path)
    parser.add_argument("--output",required=True,type=Path)
    parser.add_argument("--flpq",required=True,type=Path)
    args=parser.parse_args()
    assert not args.output.exists(),"Refusing to overwrite existing F4 test folder"
    args.output.mkdir(parents=True)
    root=args.output
    methane=args.base/"flpq_smoke"
    cr2s={x:args.base/("flpq_phase2_converged_"+x+"_12875") for x in ("on","off")}
    baseline_methane=args.base/"flpq_smoke/phase2_phase1_regression_20261016/ranks_1/LPQ"
    baselines={x:args.base/("flpq_phase2_converged_"+x+"_12875/mpi_regression/ranks_1/LPQ") for x in ("on","off")}

    # F4 default legacy must produce identical ASCII/cube with no HDF5.
    legacy=root/"default_legacy"
    prepare_case(methane,legacy)
    run_case(args.flpq,legacy,2)
    assert not any(x.endswith(".h5") for x in outputs(legacy))
    for name in outputs(legacy):
        if name.endswith((".dat",".cube")):
            assert sha256(legacy/"LPQ"/name)==sha256(baseline_methane/name)
    print("PASS default legacy SHA256 byte identical",flush=True)

    # Scoped legacy selection must not unexpectedly generate other quantities.
    legacy_selection=root/"legacy_selection"
    prepare_case(methane,legacy_selection,"Output.Quantities electron_density")
    run_case(args.flpq,legacy_selection,2)
    files=outputs(legacy_selection)
    assert files=={"LPQ.electron_density.dat","LPQ.electron_density.cube","header.cube"},files
    print("PASS legacy explicit selection",flush=True)

    # Double-precision gzip both, MPI rank variance and low-memory.
    chosen={"electron_density","charge_current"}
    combined="Output.Format both\nOutput.Quantities electron_density, charge_current\n"
    normal=root/"methane_both_n1"
    prepare_case(methane,normal,combined)
    run_case(args.flpq,normal,1)
    check_subset(normal,chosen,baseline_methane,require_legacy=True)
    multi=root/"methane_both_n4"
    prepare_case(methane,multi,combined)
    run_case(args.flpq,multi,4)
    check_subset(multi,chosen,baseline_methane,require_legacy=True)
    compare_selected(normal,multi,chosen)
    low=root/"methane_both_low_n2"
    prepare_case(methane,low,combined,low=True)
    run_case(args.flpq,low,2)
    check_subset(low,chosen,baseline_methane,require_legacy=True)
    compare_selected(normal,low,chosen)
    print("PASS methane subset MPI 1/4 and low-memory 2, legacy bytes",flush=True)

    # Atomic fail-on-existing: keep both original HDF5 files unchanged.
    names=sorted(chosen)
    hashes={n:sha256(normal/"LPQ"/("LPQ."+n+".h5")) for n in names}
    text=run_case(args.flpq,normal,1,expect_ok=False)
    assert "exists" in text.lower(),text[-800:]
    assert hashes=={n:sha256(normal/"LPQ"/("LPQ."+n+".h5")) for n in names}
    assert not list((normal/"LPQ").glob("*.partial.h5"))
    print("PASS collision preflight refuses without touching existing HDF5",flush=True)

    # Explicit replace writes new float32 file only on successful completion.
    inp=normal/"flpq.inp"
    inp.write_text(inp.read_text()+"\nOutput.HDF5.OnExisting replace\n"
                   +"Output.HDF5.Precision float32\n"
                   +"Output.HDF5.Compression none\n")
    run_case(args.flpq,normal,2)
    check_subset(normal,chosen,baseline_methane,precision="float32",
                 compression="none",require_legacy=True)
    assert any(sha256(normal/"LPQ"/("LPQ."+n+".h5"))!=hashes[n] for n in names)
    print("PASS replace atomic final swap, new float32/none, legacy intact",flush=True)

    # Orphaned partial protects against accidental resume or overwrite.
    pending=root/"stale_partial"
    prepare_case(methane,pending,"Output.Format hdf5\nOutput.Quantities electron_density")
    (pending/"LPQ").mkdir()
    partial=pending/"LPQ/LPQ.electron_density.partial.h5"
    partial.write_bytes(b"existing partial, preserve exactly")
    run_case(args.flpq,pending,2,expect_ok=False)
    assert partial.read_bytes()==b"existing partial, preserve exactly"
    assert not (pending/"LPQ/LPQ.electron_density.h5").exists()
    print("PASS stale partial rejection and preservation",flush=True)

    # Production HDF5 only, lossless no compression, gzip 1 and 9, float32.
    variants=[
        ("none64","Output.HDF5.Compression none","float64","none",4),
        ("gzip1","Output.HDF5.Compression gzip\nOutput.HDF5.GzipLevel 1","float64","gzip",1),
        ("gzip9","Output.HDF5.Compression gzip\nOutput.HDF5.GzipLevel 9","float64","gzip",9),
        ("float32","Output.HDF5.Precision float32","float32","gzip",4),
    ]
    for label,extra,precision,compression,level in variants:
        case=root/label
        prepare_case(methane,case,"Output.Format hdf5\nOutput.Quantities electron_density, charge_current\n"+extra)
        run_case(args.flpq,case,2)
        check_subset(case,chosen,baseline_methane,precision,compression,level)
    print("PASS hdf5-only precision/compression 4 scenarios",flush=True)

    for soc in ("on","off"):
        ref=baselines[soc];src=cr2s[soc]
        sel={"stress","stress_S","stress_A","charge_current"}
        config="Output.Format both\nOutput.Quantities stress,stress_S,stress_A,charge_current\n"
        full=root/("cr2_"+soc+"_n1")
        prepare_case(src,full,config)
        run_case(args.flpq,full,1)
        check_subset(full,sel,ref,require_legacy=True)
        parallel=root/("cr2_"+soc+"_n4")
        prepare_case(src,parallel,config)
        run_case(args.flpq,parallel,4)
        check_subset(parallel,sel,ref,require_legacy=True)
        compare_selected(full,parallel,sel)
        lowcr=root/("cr2_"+soc+"_low_n2")
        prepare_case(src,lowcr,config,low=True)
        run_case(args.flpq,lowcr,2)
        check_subset(lowcr,sel,ref,require_legacy=True)
        compare_selected(full,lowcr,sel)
        print("PASS Cr2 SOC",soc,"stress decomposition selection MPI/lowmem",flush=True)

    # Reject invalid selections and settings before producing science results.
    negatives=[
        ("unknown","Output.Quantities not_a_physical_quantity"),
        ("duplicate","Output.Quantities electron_density,electron_density"),
        ("invalid_spin","Output.Quantities zeta_potential"),
        ("mixed_all","Output.Quantities all,electron_density"),
        ("no_values","Output.Quantities    # comment"),
        ("duplicate_key","Output.Quantities electron_density\nOutput.Quantities charge_current"),
        ("bad_mode","Output.Format hf5"),
        ("duplicate_mode","Output.Format hdf5\nOutput.Format both"),
        ("bad_precision","Output.HDF5.Precision fp16"),
        ("bad_compression","Output.HDF5.Compression lz4"),
        ("bad_level_negative","Output.HDF5.GzipLevel -1"),
        ("bad_level_high","Output.HDF5.GzipLevel 10"),
        ("invalid_integer","Output.HDF5.GzipLevel xyz"),
        ("bad_existing","Output.HDF5.OnExisting append"),
        ("bad_chunk_y_zero","Output.HDF5.ChunkY 0"),
        ("bad_chunk_y_high","Output.HDF5.ChunkY 129"),
        ("bad_chunk_z_zero","Output.HDF5.ChunkZ 0"),
        ("bad_chunk_z_high","Output.HDF5.ChunkZ 513"),
        ("bad_block_low","Output.HDF5.BlockPoints 15"),
        ("bad_block_high","Output.HDF5.BlockPoints 8193"),
        ("bad_perf_flag","Output.HDF5.LogPerformance invalid"),
    ]
    for label,invalid in negatives:
        case=root/("invalid_"+label)
        prefix="" if label=="bad_mode" else "Output.Format hdf5\n"
        prepare_case(methane,case,prefix+invalid)
        log=run_case(args.flpq,case,2,expect_ok=False)
        assert "error" in log.lower() or "abort" in log.lower() or "quantity" in log.lower(),(label,log[-1000:])
        assert not list((case/"LPQ").glob("*.h5")) if (case/"LPQ").is_dir() else True
    print("PASS invalid configurations",len(negatives),flush=True)
    print("F4 REGRESSION PASS",flush=True)


if __name__=="__main__":
    main()
