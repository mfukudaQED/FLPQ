#!/usr/bin/env python3
"""Bounded-block edge case: BlockPoints<Nz, row batching, valid layout."""
import argparse
from pathlib import Path
from benchmark_hdf5_f5 import run_case
def main():
    p=argparse.ArgumentParser()
    p.add_argument('--fixture',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--exe',type=Path,required=True)
    p.add_argument('--baseline',type=Path,required=True)
    p.add_argument('--timeout',type=int,default=180)
    args=p.parse_args()
    assert not args.out.exists()
    args.out.mkdir(parents=True)
    case={"name":"24_block16","grid":24,"ranks":4,"block":16,
          "chunk_y":16,"chunk_z":128,"compression":"gzip","gzip_level":4}
    result=run_case(args,case,args.out/case['name'],args.baseline)
    print("F5 small-block boundary PASS:",[(x['quantity'],x['file_bytes']) for x in result])
if __name__=="__main__":main()
