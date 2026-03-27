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
                print(f'{n}\t{entry.name}')
                seq[entry.name.replace('.fa', '')] = []

    directories = [entry.name for entry in os.scandir('.') if entry.is_dir()]
    dirlist = []
    for dir in directories:
        dirlist.append(dir)
        fpt = dir + '/fpt'

        with os.scandir(fpt) as entries:
            fptcount = 0
            for entry in entries:
                if entry.name.endswith('.fpt'):
                    fptcount += 1
                    qname = entry.name.replace('.fpt', '')
                    seq[qname].append(dir)

            print(f'{dir}\t{fptcount}')

    for query in seq:
        print(f'{query}\t{len(seq[query])}')

    exit(0)
