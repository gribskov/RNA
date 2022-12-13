"""=================================================================================================
Calculate the distance between fingerprints.
Assume fingerprints are all in the same directory
Calculate the jaccard similarity and Cray-Curtis distance between all pairs

Michael Gribskov     04 February 2022
================================================================================================="""
import glob
import os
from os.path import basename
import sys
from fingerprint import Fingerprint, FingerprintSet
import datetime


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
                    default='selected.motiflist')
    args = cl.parse_args()

    # if output filename is auto, create a file name from the fingerprint file name
    if str(args.distance) == 'auto':
        out = args.distance
        out = out.replace('*', '_all')
        args.distance = out
    if not args.indir.endswith('/'):
        args.indir += '/'

    return args


def get_selected(fname):
    """---------------------------------------------------------------------------------------------
    Read a list of selected motifs from the command line (opt.motif). if fname cannot be opened for
    reading, assume the list should be all motifs. if the motif file cannot be read an empty list is
    returned. this will interpreted downstream to mean use all motifs

    :param fname: string    file of selected motif names
    :return: list           motif names
    ---------------------------------------------------------------------------------------------"""
    if not os.access(fname, os.R_OK):
        sys.stderr.write(f'{newline}motif file {fname} is not accesible{newline}')
        sys.stderr.write(f'All motifs will be used{newline}')
        return []

    # the file was openable so read it
    motifs = open(opt.motif, 'r')
    mlist = []
    for line in motifs:
        mname, count = line.rstrip().split()
        mlist.append(mname)

    motifs.close()
    return mlist


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    newline = '\n'

    opt = process_command_line()
    sys.stderr.write(f"fingerprint distance: calculate distance between sets of fingerprints{newline}")
    daytime = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S{newline}')
    sys.stderr.write(f"{newline}fingerprint files: {opt.indir}{opt.fpt}{newline}")
    sys.stderr.write(f"motif file: {opt.motif}{newline}")
    sys.stderr.write(f"distance output: {opt.distance}{newline}")

    # make a list of all fingerprints in the target directory and read into a FingerprintSet
    fpt_list = glob.glob(f'{opt.indir}{opt.fpt}')

    fpt = FingerprintSet()
    for this_fpt in fpt_list:
        f = Fingerprint()
        f.readYAML(this_fpt)
        fpt.append(f)

    fpt_n = len(fpt)
    print(f'{fpt_n} fingerprints read')

    # select only motifs in opt.motif, or use all if not provided
    motifs_selected = get_selected(opt.motif)
    fpt.select(motifs_selected)

    # distance calculation
    maximum = {'jaccard': 0, 'bray-curtis': 0}
    minimum = {'jaccard': 1000000, 'bray-curtis': 1000000}

    jaccard = fpt.jaccard_sim([])
    bc = fpt.bray_curtis_dis([])

    try:
        out = open(opt.distance, 'w')
    except OSError:
        sys.stderr.write(f'{newline} output distance file ({opt.distance}) could not be opened{newline}')
        exit(1)

    for d in jaccard:
        i = d[0]
        j = d[1]
        out.write(f'{basename(fpt[i].information["File"])}\t')
        out.write(f'{basename(fpt[j].information["File"])}\t')
        out.write(f'{d[2]:.3f}\t')
        out.write(f'{bc[j][2]:.3f}\n')

        maximum['jaccard'] = max(maximum['jaccard'], d[2])
        minimum['jaccard'] = min(minimum['jaccard'], d[2])
        maximum['bray-curtis'] = max(maximum['bray-curtis'], bc[j][2])
        minimum['bray-curtis'] = min(minimum['bray-curtis'], bc[j][2])

    print()
    for dist in ('jaccard', 'bray-curtis'):
        print(f'{dist}\tmaximum={maximum[dist]:.3f}\tminimum={minimum[dist]:.3f}')

    exit(0)
