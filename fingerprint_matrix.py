"""=================================================================================================
Convert a set of fingerprint files into a binary form for rapid reading and manipultation
================================================================================================="""
from fingerprint import FingerprintMatrix

# main program
if __name__ == '__main__':
    selection = 'data/fpt/*'
    fmat = FingerprintMatrix()
    fmat.read_files(selection)

    exit(0)
