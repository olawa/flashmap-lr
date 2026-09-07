#!/usr/bin/env python3
"""Generate independent synthetic train/heldout exact 5kb reads for MAPQ plumbing.
Not a sequencing-error model or biological calibration dataset.
Usage: python3 scripts/mapq_fixture.py OUTPUT_DIRECTORY
"""
import pathlib
import random
import sys

root = pathlib.Path(sys.argv[1])
root.mkdir(parents=True, exist_ok=True)
for split, seed in [('train', 719), ('heldout', 991)]:
    rng = random.Random(seed)
    repeated = ''.join(rng.choices('ACGT', k=25000))
    unique = ''.join(rng.choices('ACGT', k=25000))
    (root / f'{split}.fa').write_text(f'>a\n{repeated}\n>b\n{repeated}\n>unique\n{unique}\n')
    reads, truth = [], []
    for i in range(600):
        contig = ['a', 'b', 'unique'][i % 3]
        start = rng.randrange(100, 19000)
        sequence = (unique if contig == 'unique' else repeated)[start:start+5000]
        name = f'{split}_{i}'
        reads.append(f'>{name}\n{sequence}\n')
        truth.append(f'{name}\t{contig}\t{start}\t{start+5000}\t+\n')
    (root / f'{split}.reads.fa').write_text(''.join(reads))
    (root / f'{split}.truth.tsv').write_text(''.join(truth))
