"""=================================================================================================
Convert a set of fingerprint files into a binary form for rapid reading and manipulation
================================================================================================="""
from fingerprint import FingerprintMatrix

# main program
if __name__ == '__main__':
    selection = 'data/fpt/*'
    fmat = FingerprintMatrix()
    fmat.read_files(selection)
    fmat.select_min_max(2, 2)
    fmat.select_min_max(0, 7, setval=True, recalculate=True)
    # fmat.index2matrix()

    exit(0)
