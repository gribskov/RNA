"""=================================================================================================
Calculate the distance between fingerprints.
Assume fingerprints are all in the same directory
Calculate the jaccard similarity and Cray-Curtis distance between all pairs

Michael Gribskov     04 February 2022
================================================================================================="""
import glob
from os.path import basename
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
    cl.add_argument('-m', '--indir',
                    help='Directory with fingerprint files (default=%(default)s)',
                    default='../data/fpt')
    cl.add_argument('-r', '--fpt ',
                    help='filename for fingerprint files, can be wildcard  (default=%(default)s',
                    default='.fpt')
    cl.add_argument('-o', '--distance',
                    help='Distance output file name (default=%(default)s, auto=calculate from fpt)',
                    default='fingerprint.dist')
    cl.add_argument('-m', '--outdir',
                    help='Directory for output distance files (default=%(default)s)',
                    default='./')
    args = cl.parse_args()

    # if output filename is auto, create a file name from the fingerprint file name
    if str(args.distance) == 'auto':
        out = args.distance
        out = out.replace('*', '_all')
        args.distance = out
    if args.dir.endswith('/'):
        pass
    else:
        args.dir += '/'

    return args


# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    daytime = datetime.datetime.now()

    # setup directories and suffix according to command line
    opt = process_command_line()
    fpt_dir = opt.fpt
    fpt_suffix = opt.suffix
    outfile = opt.out

    # make a list of all fingerprints in the target directory
    fpt_list = glob.glob(f'{fpt_dir}*{fpt_suffix}')

    fpt = FingerprintSet()
    for this_fpt in fpt_list:
        f = Fingerprint()
        f.readYAML(this_fpt)
        fpt.append(f)

    fpt_n = len(fpt)
    print(f'{fpt_n} fingerprints read')

    # distance calculation
    maximum = {'jaccard': 0, 'bray-curtis': 0}
    minimum = {'jaccard': 1000000, 'bray-curtis': 1000000}

    jaccard = fpt.jaccard_sim([])
    bc = fpt.bray_curtis_dis([])

    out = open(outfile, 'w')
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
