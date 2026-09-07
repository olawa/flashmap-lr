#!/usr/bin/env python3
"""Deterministic 1000-read end-to-end comparison; verifies SAM record identity."""
import random, pathlib, subprocess, time, hashlib, json, statistics
import argparse
p=argparse.ArgumentParser();p.add_argument('--baseline',required=True);p.add_argument('--optimized',required=True);p.add_argument('--out',required=True);args=p.parse_args()
root=pathlib.Path(args.out);root.mkdir(parents=True,exist_ok=True);rng=random.Random(906)
ref=bytearray(rng.choices(b'ACGT',k=400000))
for start in range(1500,len(ref)-150,2500): ref[start:start+120]=b'ACGTAC'*20
(root/'ref.fa').write_bytes(b'>ref\n'+ref+b'\n')
with (root/'reads.fa').open('wb') as out:
 for i in range(1000):
  start=rng.randrange(len(ref)-12000);read=bytearray(ref[start:start+10000]);j=50
  while j<len(read)-50:
   if rng.random()<0.007:
    n=rng.choice([1,2,3,4,8,20,40])
    if rng.random()<0.5: del read[j:j+n]
    else: read[j:j]=bytes(rng.choices(b'ACGT',k=n));j+=n
   elif rng.random()<0.002:read[j]=rng.choice(b'ACGT')
   j+=1
  if i%2:read=read.translate(bytes.maketrans(b'ACGT',b'TGCA'))[::-1]
  out.write(f'>read{i}\n'.encode()+read+b'\n')
rows=[]
for workers in [1,4]:
 for rep in range(3):
  for label,binary in [('corrected',args.baseline),('optimized',args.optimized)]:
   sam=root/f'{label}-{workers}.sam'; log=root/f'{label}-{workers}.log'
   cmd=[binary,'--reference',str(root/'ref.fa'),'--reads',str(root/'reads.fa'),'--output',str(sam),'--workers',str(workers),'--dual-affine','--profile']
   t=time.perf_counter()
   result=subprocess.run(cmd,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE,check=True)
   elapsed=time.perf_counter()-t;log.write_bytes(result.stderr)
   records=b''.join(line for line in sam.read_bytes().splitlines(keepends=True) if not line.startswith(b'@'))
   row=dict(label=label,workers=workers,rep=rep,seconds=elapsed,sha256=hashlib.sha256(records).hexdigest());rows.append(row);print(json.dumps(row),flush=True)
(root/'results.json').write_text(json.dumps(rows,indent=2))

assert len({row["sha256"] for row in rows}) == 1, "alignment output changed"
for workers in [1,4]:
 med={label:statistics.median(r["seconds"] for r in rows if r["workers"]==workers and r["label"]==label) for label in ["corrected","optimized"]}
 print(workers, med, "speedup", med["corrected"]/med["optimized"])
