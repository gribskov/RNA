"""=================================================================================================
Calculate the distance between fingerprints.
Assume fingerprints are all in the same directory
Calculate the jaccard similarity and Bray-Curtis distance between all pairs

Michael Gribskov     04 February 2022
================================================================================================="""
from datetime import datetime
import glob
import sys
from os.path import basename

from fingerprint import Fingerprint, FingerprintSet


def process_command_line():
    """---------------------------------------------------------------------------------------------
    read command line options and return as a Namespace object. The Namespace object behaves much as
    a dictionary, and can be converted to a dictionary using the vars() method

    :return: Namespace object
    ---------------------------------------------------------------------------------------------"""
    import argparse

    cl = argparse.ArgumentParser(
        description='Calculate distances between fingerprints',
        formatter_class=lambda prog: argparse.HelpFormatter(prog, width=120, max_help_position=40)
        )
    cl.add_argument('-i', '--indir',
                    help='Directory with fingerprint files (default=%(default)s)',
                    default='../data/fpt/')
    cl.add_argument('-f', '--fpt',
                    help='filename for fingerprint files, can be wildcard  (default=%(default)s)',
                    default='*.fpt')
    cl.add_argument('-d', '--distance',
                    help='Distance output file name (default=%(default)s, auto=calculate from fpt)',
                    default='fingerprint.distance')
    cl.add_argument('-o', '--outdir',
                    help='Directory for output distance files (default=%(default)s)',
                    default='./')
    cl.add_argument('-m', '--motif',
                    help='Selected motif file (default=%(default)s)',
                    default='')
    args = cl.parse_args()

    # if output filename is auto, create a file name from the fingerprint file name
    if str(args.distance) == 'auto':
        out = args.distance
        out = out.replace('*', '_all')
        args.distance = out
    if not args.indir.endswith('/'):
        args.indir += '/'

    return args


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    opt = process_command_line()
    now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    sys.stderr.write(f"fingerprint distance: calculate distance between sets of fingerprints\t {now}\n")
    sys.stderr.write(f'')
    sys.stderr.write(f"\nfingerprint files: {opt.indir}{opt.fpt}\n")
    sys.stderr.write(f"distance output: {opt.distance}\n")
    if opt.motif:
        sys.stderr.write(f"motif file: {opt.motif}\n")
    else:
        sys.stderr.write(f"Use all motifs\n")

    try:
        out = open(opt.distance, 'w')
    except OSError:
        sys.stderr.write(f'\n output distance file ({opt.distance}) could not be opened\n')
        exit(1)

    # make a list of all fingerprints in the target directory and read into a FingerprintSet
    fpt_list = glob.glob(f'{opt.indir}{opt.fpt}')

    fpt = FingerprintSet()
    for this_fpt in fpt_list:
        f = Fingerprint()
        f.readYAML(this_fpt)
        fpt.append(f)

    fpt_n = len(fpt)
    print(f'\n{fpt_n} fingerprints read from {opt.indir}{opt.fpt}')

    # select only motifs in opt.motif, or use all if not provided and construct binary matrix
    # set binary motif matrix
    motif_n = fpt.select(filename=opt.motif)
    fpt.binary_matrix()
    print(f'{len(fpt.i2motif)} motifs selected from {motif_n}')

    # distance calculation
    # jaccard = fpt.jaccard_binary()
    jaccard = fpt.jaccard_scale()
    bc = fpt.bray_curtis_binary()

    for d in jaccard:
        i = d[0]
        j = d[1]
        out.write(f'{basename(fpt[i].information["File"])}\t')
        out.write(f'{basename(fpt[j].information["File"])}\t')
        out.write(f'{d[2]:.3f}\t')
        out.write(f'{bc[j][2]:.3f}\n')

    jmax = max(jaccard, key=lambda a: a[2])
    jmin = min(jaccard, key=lambda a: a[2])
    print(f'\nJaccard similarity:\tmaximum={jmax[2]:.3f}\tminimum={jmin[2]:.3f}')
    bmax = max(bc, key=lambda a: a[2])
    bmin = min(bc, key=lambda a: a[2])
    print(f'Bray-Curtis distance:\tmaximum={bmax[2]:.3f}\tminimum={bmin[2]:.3f}')

    exit(0)
