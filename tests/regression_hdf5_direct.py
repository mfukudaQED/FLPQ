#!/usr/bin/env python3
"""F3 direct-HDF5 acceptance tests on already-computed OpenMX fixtures.

Run in an allocated compute node with a functioning mpiexec launcher.
Never touches the fixture directories. Requires numpy/h5py.
"""
import argparse
import hashlib
from pathlib import Path
import shutil
import subprocess
import os
import sys

import h5py
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"tools"))
from flpq_convert import FLPQFile

TENSOR_ORDER=(0,3,6,1,4,7,2,5,8)


def sha256(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for piece in iter(lambda: f.read(1<<20),b""):
            h.update(piece)
    return h.hexdigest()


def rows(path):
    return np.asarray([
        [float(x) for x in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ],dtype=float)


def make_case(fixture,out, mode,low):
    if out.exists():
        raise AssertionError("Refusing to overwrite "+str(out))
    out.mkdir(parents=True)
    shutil.copy2(fixture/"OpenMX_LPQ.bin",out/"OpenMX_LPQ.bin")
    shutil.copytree(fixture/"OpenMX_LPQ_DM",out/"OpenMX_LPQ_DM")
    txt=(fixture/"flpq.inp").read_text()
    assert "Output.Format" not in txt
    if low:
        assert "save_memory_3d off" in txt
        txt=txt.replace("save_memory_3d off","save_memory_3d on")
    (out/"flpq.inp").write_text(txt+"\nOutput.Format "+mode+"\n")


def run(exe, folder,rank,env):
    cmd=["mpiexec","-n",str(rank),str(exe),"flpq.inp"]
    with (folder/"flpq.log").open("w") as log:
        p=subprocess.run(cmd,cwd=str(folder),env=env,stdout=log,stderr=subprocess.STDOUT,timeout=270)
    if p.returncode:
        print((folder/"flpq.log").read_text()[-4000:],flush=True)
        raise AssertionError("FLPQ failed with exit "+str(p.returncode)+" in "+str(folder))


def field_array(h5):
    if h5.quantity=="stress_diag":
        ev=h5.file["/field/eigenvalues"][:]
        v=h5.file["/field/eigenvectors"][:]
        vectors=np.moveaxis(v,(0,1),(-2,-1)).swapaxes(-2,-1).reshape(-1,9)
        return np.concatenate((ev.reshape(3,-1).T,vectors),axis=1)
    a=h5.file["/field/values"][:]
    result=a.reshape(h5.ncomp,-1).T
    if h5.field_rank==2: return result[:,TENSOR_ORDER]
    return result


def check_case(folder, legacy_reference=None):
    lpq=folder/"LPQ"
    h5s=sorted(lpq.glob("LPQ.*.h5"))
    assert len(h5s)>=9,(folder,len(h5s))
    dat={p.name.removesuffix(".dat").removeprefix("LPQ."):p
         for p in lpq.glob("*.dat")}
    for hp in h5s:
        with FLPQFile(hp) as h:
            assert len(h.shape)==3 and h.natoms>0
            assert np.all(np.isfinite(h.origin)) and np.all(np.isfinite(h.steps))
            data=field_array(h)
            name=h.quantity
            if dat:
                assert name in dat
                original=rows(dat[name])
                assert original.shape==(np.prod(h.shape), 3+h.ncomp),(name,original.shape,h.ncomp)
                np.testing.assert_allclose(data,original[:,3:],rtol=6e-9,atol=1e-12,
                                           err_msg="HDF5 vs legacy "+name)
                for ix in (0,original.shape[0]//2,original.shape[0]-1):
                    i,j,k=np.unravel_index(ix,h.shape)
                    xyz=h.origin+i*h.steps[0]+j*h.steps[1]+k*h.steps[2]
                    np.testing.assert_allclose(original[ix,:3],xyz,rtol=8e-5,atol=2e-4,
                                               err_msg="grid geometry "+name)
    if dat:
        assert len(dat)==len(h5s)
    else:
        assert not list(lpq.glob("*.cube"))
    if legacy_reference:
        old={p.name:p for p in legacy_reference.iterdir()
             if p.suffix in (".dat",".cube")}
        current={p.name:p for p in lpq.iterdir()
                 if p.suffix in (".dat",".cube")}
        assert current.keys()==old.keys(),(folder,current.keys(),old.keys())
        for name,p in old.items():
            assert sha256(current[name])==sha256(p),"Legacy byte mismatch "+name
    return len(h5s)


def compare_hdf5(a,b):
    for p in sorted((a/"LPQ").glob("LPQ.*.h5")):
        other=b/"LPQ"/p.name
        assert other.exists(),other
        with FLPQFile(p) as x, FLPQFile(other) as y:
            assert x.shape==y.shape
            assert x.quantity==y.quantity
            np.testing.assert_allclose(field_array(x),field_array(y),
                                       rtol=3e-13,atol=1e-13,
                                       err_msg="MPI/low-memory "+p.name)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--flpq",type=Path,required=True)
    parser.add_argument("--fixture-base",type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():parser.error("Output directory must not exist")
    args.output.mkdir(parents=True)
    env=os.environ.copy()
    results={}
    sources=[
        ("methane",args.fixture_base/"flpq_smoke",args.fixture_base/"flpq_smoke/phase2_phase1_regression_20261016/ranks_1/LPQ"),
        ("cr2_on",args.fixture_base/"flpq_phase2_converged_on_12875",args.fixture_base/"flpq_phase2_converged_on_12875/mpi_regression/ranks_1/LPQ"),
        ("cr2_off",args.fixture_base/"flpq_phase2_converged_off_12875",args.fixture_base/"flpq_phase2_converged_off_12875/mpi_regression/ranks_1/LPQ"),
    ]
    for label,fixture,baseline in sources:
        maincase=None
        for n in (1,2,4):
            case=args.output/(label+"_n"+str(n))
            make_case(fixture,case,"both",False)
            run(args.flpq,case,n,env)
            nfields=check_case(case,baseline)
            if maincase:compare_hdf5(maincase,case)
            else:maincase=case
            print("PASS",label,"MPI",n,"fields",nfields,flush=True)
        low=args.output/(label+"_lowmem_n2")
        make_case(fixture,low,"both",True)
        run(args.flpq,low,2,env)
        check_case(low)
        compare_hdf5(maincase,low)
        # low-memory legacy dat should match full-grid legacy dat.
        for p in (maincase/"LPQ").glob("*.dat"):
            assert sha256(p)==sha256(low/"LPQ"/p.name),"Low memory dat differs "+p.name
        print("PASS",label,"lowmem n2",flush=True)
        results[label]=nfields
    # HDF5-only mode: no ASCII.
    fixture=args.fixture_base/"flpq_smoke"
    pure=args.output/"methane_hdf5_only"
    make_case(fixture,pure,"hdf5",False)
    run(args.flpq,pure,2,env)
    check_case(pure)
    compare_hdf5(args.output/"methane_n1",pure)
    print("PASS methane HDF5-only no legacy artifacts",flush=True)
    print("F3 PASS:",results,"HDF5-only, MPI 1/2/4, low memory, legacy SHA256",flush=True)


if __name__=="__main__":
    main()
