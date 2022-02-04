"""=================================================================================================
Calculate the distance between fingerprints.
Assume fingerprint are all in the same directory
Calculate the jaccard distance between all pairs

Michael Gribskov     04 February 2022
================================================================================================="""
import glob
from fingerprint import Fingerprint, FingerprintSet

# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    fpt_dir = './'
    fpt_suffix = '.fpt'

    # make a list of all fingerprints in the target directory
    fpt_list = glob.glob(f'{fpt_dir}*{fpt_suffix}')

    fpt = FingerprintSet()
    for this_fpt in fpt_list:
        f = Fingerprint()
        f.readYAML(this_fpt)
        fpt.append(f)


    # distance calculation
    exit(0)
