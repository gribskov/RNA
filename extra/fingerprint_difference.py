"""=====================================================================================================================
fingerprint_difference.py
compare two fingerprints and report on the difference in motifs

2026-04-14 gribskov
====================================================================================================================="""
import sys

from fingerprint import Fingerprint

# ======================================================================================================================
# Main
# ======================================================================================================================
if __name__ == '__main__':
    fpt1 = open(sys.argv[1], 'r')
    fpt2 = open(sys.argv[2], 'r')

    f1 = Fingerprint()
    f1.readYAML(fpt1)
    m = {k:1 for k in f1.motif}

    f2 = Fingerprint()
    f2.readYAML(fpt2)
    for m2 in f2.motif:
        if m2 in m:
            m[m2] += 2
        else:
            m[m2] =2

    n1 = 0
    n2 = 0
    nboth = 0
    for motif in m:
        if m[motif] == 1:
            n1 += 1
        elif m[motif] == 2:
            n2 += 1
        else:
            n1 += 1
            n2 += 1
            nboth += 1

    print(f'motifs in {sys.argv[1]}: {n1}')
    print(f'motifs in {sys.argv[2]}: {n2}')
    print(f'motifs in both: {nboth}')

    print(f'\nmotifs only in {sys.argv[1]}: {n1-nboth}')
    n = 0
    for motif in sorted(m):
        if m[motif] == 1:
            n += 1
            print(f'\t{n}\t{motif}')

    print(f'\nmotifs only in {sys.argv[2]}: {n2-nboth}')
    n = 0
    for motif in sorted(m):
        if m[motif] == 2:
            n += 1
            print(f'\t{n}\t{motif}')

    exit(0)
