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
    indata = '../../RNA/data/fasta_fixed/*.fa'
    with os.scandir(indata) as entries:
        n = 0
        for entry in entries:
            if entry.name.endswith('fa'):
                n += 1
                print(f'{n}\t{entry.name}')



    directories = [entry.name for entry in os.scandir('.') if entry.is_dir()]

    for dir in directories:
        print(dir)
        fpt = dir + '/fpt'

        with os.scandir('.') as entries:
            for entry in entries:
                if entry.is_file():
                    print(f"{entry.name}: {entry.stat().st_size} bytes")

    exit(0)
