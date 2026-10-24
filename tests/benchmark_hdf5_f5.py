#!/usr/bin/env python3
"""F5 real-FLPQ writer bench on a sham compute node.

No synthetic values: reuses converged Cr2 SOC OpenMX density matrix.
Benchmarks density(1), charge current(3), stress(9).
Each case has a fresh output path; no user/global temporary locations.
"""
import argparse,csv,json,re,shutil,subprocess,sys,time
from pathlib import Path
import h5py,numpy as np

CASE_SCHEMA = ["name","grid","ranks","block","chunk_y","chunk_z","compression","gzip_level"]
PERF = re.compile(r"FLPQ_HDF5_PERF\s+(.*)$")
LOG_KEYS = ["quantity","slice","ranks","nx","ny","nz","block_points","chunk_y","chunk_z",
"precision","gzip","gzip_level","io_seconds","wall_seconds","eval_max_seconds","max_rss_kib","file_bytes"]

def setup(fixture,target,case):
    assert not target.exists(),target
    target.mkdir(parents=True)
    shutil.copy2(fixture/"OpenMX_LPQ.bin",target/"OpenMX_LPQ.bin")
    shutil.copytree(fixture/"OpenMX_LPQ_DM",target/"OpenMX_LPQ_DM")
    raw=(fixture/"flpq.inp").read_text()
    assert "Output.Format" not in raw and "Output.HDF5." not in raw
    raw,n=re.subn(r"(?im)^Dens_Ngrid\s+[0-9]+\s+[0-9]+\s+[0-9]+",
                   "Dens_Ngrid "+str(case["grid"])+" "+str(case["grid"])+" "+str(case["grid"]),raw,count=1)
    assert n==1
    raw += "\nOutput.Format hdf5\nOutput.Quantities electron_density,charge_current,stress\n"
    raw += "Output.HDF5.LogPerformance on\n"
    raw += "Output.HDF5.BlockPoints {}\n".format(case["block"])
    raw += "Output.HDF5.ChunkY {}\n".format(case["chunk_y"])
    raw += "Output.HDF5.ChunkZ {}\n".format(case["chunk_z"])
    raw += "Output.HDF5.Compression {}\n".format(case["compression"])
    raw += "Output.HDF5.GzipLevel {}\n".format(case["gzip_level"])
    (target/"flpq.inp").write_text(raw)
def parse_log(s):
    records=[]
    for line in s.splitlines():
        m=PERF.search(line)
        if m:
            d=dict(re.findall(r"([A-Za-z_]+)=([^\s]+)",m.group(1)))
            if set(LOG_KEYS)-set(d):raise ValueError("missing columns "+str(d))
            records.append(d)
    return records

def compare_outputs(a,b,precision="float64"):
    names=["electron_density","charge_current","stress"]
    for q in names:
        pa=a/"LPQ"/("LPQ."+q+".h5")
        pb=b/"LPQ"/("LPQ."+q+".h5")
        with h5py.File(str(pa),"r") as fa,h5py.File(str(pb),"r") as fb:
            x=fa["field/values"][:];y=fb["field/values"][:]
            assert x.dtype==y.dtype
            np.testing.assert_allclose(x,y,rtol=2e-12,atol=5e-13,err_msg=q)
            assert np.array_equal(fa["grid/shape"],fb["grid/shape"])
            assert np.array_equal(fa["grid/step_vectors"],fb["grid/step_vectors"])

def run_case(args,case,path,baseline):
    setup(args.fixture,path,case)
    cmd=["mpiexec","-n",str(case["ranks"]),str(args.exe),"flpq.inp"]
    start=time.monotonic()
    with (path/"flpq.log").open("w") as fd:
        r=subprocess.run(cmd,cwd=str(path),stdout=fd,stderr=subprocess.STDOUT,timeout=args.timeout)
    elapsed=time.monotonic()-start
    content=(path/"flpq.log").read_text(errors="replace")
    if r.returncode:
        raise RuntimeError("Failed {}\n{}".format(case["name"],content[-7000:]))
    rec=parse_log(content)
    assert len(rec)==3,(case["name"],len(rec),content[-4500:])
    assert {item["quantity"] for item in rec}=={"electron_density","charge_current","stress"}
    assert not list((path/"LPQ").glob("*.partial.h5"))
    assert len(list((path/"LPQ").glob("LPQ.*.h5")))==3
    if baseline is not None:compare_outputs(path,baseline)
    results=[]
    for item in rec:
        row=dict(case)
        row.update(item)
        row["name"]=case["name"]
        row["elapsed_case_seconds"]=round(elapsed,6)
        results.append(row)
    return results

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--fixture",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--exe",type=Path,required=True)
    ap.add_argument("--stage",choices=("baseline","optimized"),required=True)
    ap.add_argument("--timeout",type=int,default=300)
    args=ap.parse_args()
    if args.out.exists():ap.error("Refusing to reuse benchmark root")
    args.out.mkdir(parents=True)
    cases=[
        {"name":"24_default_n1","grid":24,"ranks":1,"block":512,"chunk_y":8,"chunk_z":32,"compression":"gzip","gzip_level":4},
        {"name":"24_default_n2","grid":24,"ranks":2,"block":512,"chunk_y":8,"chunk_z":32,"compression":"gzip","gzip_level":4},
        {"name":"24_default_n4","grid":24,"ranks":4,"block":512,"chunk_y":8,"chunk_z":32,"compression":"gzip","gzip_level":4},
        {"name":"24_gzip_none","grid":24,"ranks":4,"block":512,"chunk_y":8,"chunk_z":32,"compression":"none","gzip_level":4},
        {"name":"24_gzip_1","grid":24,"ranks":4,"block":512,"chunk_y":8,"chunk_z":32,"compression":"gzip","gzip_level":1},
        {"name":"24_gzip_9","grid":24,"ranks":4,"block":512,"chunk_y":8,"chunk_z":32,"compression":"gzip","gzip_level":9},
        {"name":"24_chunk_1x24","grid":24,"ranks":4,"block":512,"chunk_y":1,"chunk_z":24,"compression":"gzip","gzip_level":4},
        {"name":"24_chunk_24x24","grid":24,"ranks":4,"block":512,"chunk_y":24,"chunk_z":24,"compression":"gzip","gzip_level":4},
        {"name":"24_block_128","grid":24,"ranks":4,"block":128,"chunk_y":8,"chunk_z":32,"compression":"gzip","gzip_level":4},
        {"name":"24_block_2048","grid":24,"ranks":4,"block":2048,"chunk_y":8,"chunk_z":32,"compression":"gzip","gzip_level":4},
        {"name":"40_default_n4","grid":40,"ranks":4,"block":512,"chunk_y":8,"chunk_z":32,"compression":"gzip","gzip_level":4},
        {"name":"40_default_n8","grid":40,"ranks":8,"block":512,"chunk_y":8,"chunk_z":32,"compression":"gzip","gzip_level":4},
    ]
    # Compare matching-grid runs only, preserves source data and catches MPI-induced differences.
    baselines={}
    fields=list(CASE_SCHEMA)+["quantity","slice","nx","ny","nz","block_points","precision","gzip","io_seconds","wall_seconds","eval_max_seconds","max_rss_kib","file_bytes","elapsed_case_seconds"]
    csvout=args.out/"results.csv"
    with csvout.open("w",newline="") as fd:
        writer=csv.DictWriter(fd,fieldnames=fields,extrasaction="ignore")
        writer.writeheader()
        for index,case in enumerate(cases,1):
            path=args.out/case["name"]
            baseline=baselines.get(case["grid"])
            print("RUN {}/{} {} grid={} ranks={} {}".format(index,len(cases),case["name"],case["grid"],case["ranks"],case["compression"]),flush=True)
            result=run_case(args,case,path,baseline)
            writer.writerows(result)
            fd.flush()
            if baseline is None:baselines[case["grid"]]=path
            print("PASS {} sec={:.3f} sizes={} io={}".format(case["name"],float(result[0]["elapsed_case_seconds"]),
              [int(x["file_bytes"]) for x in result], [round(float(x["io_seconds"]),3) for x in result]),flush=True)
    print("F5 {} SUCCESS: {} measured configurations x 3 quantities".format(args.stage,len(cases)),flush=True)

if __name__=="__main__":main()
