"""=====================================================================================================================
check_files.py

chheck parameter sweep outputs to see which jobs never ran

Michael Gribskov 3/27/2026
====================================================================================================================="""

import os

# ======================================================================================================================
# main
# ======================================================================================================================
if __name__ == '__main__':
    # get a list of input files
    indata = '../../RNA/data/fasta_fixed'
    seq = {}
    with os.scandir(indata) as entries:
        n = 0
        for entry in entries:
            if entry.name.endswith('fa'):
                n += 1
                # print(f'{n}\t{entry.name}')
                seq[entry.name.replace('.fa', '')] = []

    print(f'{len(seq)} sequences in {indata}')

    directories = [entry.name for entry in os.scandir('.') if entry.is_dir()]
    dirlist = []
    fptlist = {}
    ecount = 0
    for dir in directories:
        dirlist.append(dir)
        fpt = dir + '/fpt'

        ecount += 1
        with os.scandir(fpt) as entries:
            fptcount = 0
            for entry in entries:
                if entry.name.endswith('.fpt'):
                    fptcount += 1
                    qname = entry.name.replace('.fpt', '')
                    seq[qname].append(dir)

            fptlist[dir] = fptcount

    # list of count of fingerprints for each parameter set by count
    for d in sorted(fptlist, key = lambda d: fptlist[d], reverse = True):
        print(f'{d}\t{fptlist[d]}')

    print(f'{ecount} parameter conditions')
    for query in sorted(seq, key=lambda q: len(seq[q]), reverse=True):
        print(f'{query}\t{len(seq[query])}')

    # missing
    for query in sorted(seq, key=lambda q: len(seq[q]), reverse=True):
        print(f'{query}\t{len(seq[query])}')
        for d in dirlist:
            if d in seq[query]:
                pass
                # print(f'{d} is in list')
            else:
                print(f'\t{d}')

    exit(0)
