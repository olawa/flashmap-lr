#!/usr/bin/env python3
"""Fit conservative MAPQ caps or evaluate primary SAM against source-locus truth.
Truth: tab-separated read, contig, start, end, strand (+/-), without header;
coordinates are zero-based, half-open. One source locus per read is required.
Use independent reference families for fitting and evaluation. SAM must contain
all reads, including unmapped records; secondary/supplementary are ignored.
"""
import argparse
import hashlib
import json
import math
import re
from pathlib import Path


def upper_error(errors, count):
    """Upper endpoint of the two-sided 95% Wilson interval."""
    if not count:
        return 1.0
    z = 1.959963984540054
    p = errors / count
    return (p + z*z/(2*count) + z*math.sqrt(p*(1-p)/count + z*z/(4*count*count))) / (1+z*z/count)


def observations(sam, truth, tolerance):
    loci = {}
    for line in Path(truth).read_text().splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        name, contig, start, end, strand = line.split('\t')
        if name in loci or strand not in ('+', '-') or int(start) < 0 or int(end) <= int(start):
            raise ValueError('invalid or duplicate truth locus: ' + name)
        loci[name] = (contig, int(start), int(end), strand)
    counts, errors = [0]*61, [0]*61
    seen, unmapped = set(), 0
    with open(sam) as stream:
        for line in stream:
            if line.startswith('@'):
                continue
            f = line.rstrip().split('\t')
            flag = int(f[1])
            if flag & 0x900:
                continue
            name = f[0]
            if name not in loci or name in seen:
                raise ValueError('unknown or duplicate primary read: ' + name)
            seen.add(name)
            if flag & 4:
                unmapped += 1
                continue
            q = int(f[4])
            if not 0 <= q <= 60:
                raise ValueError('expected rs-lra MAPQ 0..60')
            ops = re.findall(r'([1-9][0-9]*)([MIDNSHP=X])', f[5])
            if ''.join(n+op for n, op in ops) != f[5] or not ops:
                raise ValueError('invalid mapped CIGAR')
            start = int(f[3])-1
            end = start + sum(int(n) for n, op in ops if op in 'MDN=X')
            contig, ts, te, strand = loci[name]
            correct = f[2] == contig and ('-' if flag & 16 else '+') == strand and abs(start-ts) <= tolerance and abs(end-te) <= tolerance
            counts[q] += 1
            errors[q] += not correct
    if seen != loci.keys():
        raise ValueError(f'{len(loci.keys()-seen)} truth reads missing from SAM')
    if not sum(counts):
        raise ValueError('no mapped reads to calibrate/evaluate')
    return counts, errors, unmapped


def fit(counts, errors):
    # Unobserved scores get zero support. Propagating suffix minima downward
    # enforces monotonicity without increasing any statistically estimated cap.
    caps = [min(q, max(0, math.floor(-10*math.log10(upper_error(e, n)))))
            for q, (n, e) in enumerate(zip(counts, errors))]
    for q in range(59, -1, -1):
        caps[q] = min(caps[q], caps[q+1])
    return caps


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['fit', 'evaluate'])
    parser.add_argument('--sam', required=True)
    parser.add_argument('--truth', required=True)
    parser.add_argument('--tolerance', type=int, default=100)
    parser.add_argument('--table', required=True, help='output for fit; input for evaluate')
    parser.add_argument('--label', default='unspecified dataset; not transferable without validation')
    args = parser.parse_args()
    if args.tolerance < 0:
        parser.error('tolerance must be nonnegative')
    counts, errors, unmapped = observations(args.sam, args.truth, args.tolerance)
    if args.command == 'fit':
        caps = fit(counts, errors)
        digest = hashlib.sha256(Path(args.truth).read_bytes()).hexdigest()
        Path(args.table).write_text(f'# {args.label}\n# truth_sha256={digest} tolerance={args.tolerance}\n# 95% Wilson upper error caps; validate on independent loci\n' + ''.join(f'{q}\t{cap}\n' for q, cap in enumerate(caps)))
    else:
        rows = [line.split() for line in Path(args.table).read_text().splitlines() if line.strip() and not line.startswith('#')]
        if len(rows) != 61 or any(len(row) != 2 or int(row[0]) != q for q, row in enumerate(rows)):
            raise ValueError('table must contain ordered rows 0..60')
        caps = [int(row[1]) for row in rows]
        if any(not 0 <= cap <= q or (q and cap < caps[q-1]) for q, cap in enumerate(caps)):
            raise ValueError('invalid calibration caps')
    print(json.dumps({'mapped': sum(counts), 'unmapped': unmapped, 'errors': sum(errors),
        'raw_expected_errors': sum(n*10**(-q/10) for q, n in enumerate(counts)),
        'calibrated_expected_errors': sum(n*10**(-caps[q]/10) for q, n in enumerate(counts)),
        'bins': [{'raw_mapq': q, 'cap': caps[q], 'count': n, 'errors': errors[q], 'error_upper_95': upper_error(errors[q], n)} for q, n in enumerate(counts) if n]}, indent=2))

if __name__ == '__main__':
    main()
