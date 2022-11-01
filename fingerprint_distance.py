"""=================================================================================================
Calculate the distance between fingerprints.
Assume fingerprint are all in the same directory
Calculate the jaccard distance between all pairs

Michael Gribskov     04 February 2022
================================================================================================="""
import argparse
import datetime
import glob
from os.path import basename

from fingerprint import Fingerprint, FingerprintSet


def process_command_line():
    """---------------------------------------------------------------------------------------------
    command line arguments using argparse
    ---------------------------------------------------------------------------------------------"""
    cl = argparse.ArgumentParser(
        description='Calculate distances between fingerprints',
        formatter_class=lambda prog: argparse.HelpFormatter(prog, width=120, max_help_position=40)
    )
    cl.add_argument('-d', '--dir',
                    help='Directory with fingerprint files (default=%(default)s)',
                    default='../data/fpt')
    cl.add_argument('-s', '--suffix',
                    help='Suffix for fingerprint files  (default=%(default)s)',
                    default='.fpt')
    cl.add_argument('-o', '--out',
                    help='Output file name (default=%(default)s)',
                    default='distance.out')

    args = cl.parse_args()
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

    opt = process_command_line()
    outfile = opt.out
    fpt_dir = opt.dir
    fpt_suffix = opt.suffix

    # make a list of all fingerprints in the target directory
    fpt_list = glob.glob(f'{fpt_dir}*{fpt_suffix}')

    fpt = FingerprintSet()
    fpt_n = 0
    for this_fpt in fpt_list:
        fpt_n += 1
        print(f'{fpt_n}\t processing{this_fpt}')
        f = Fingerprint()
        f.readYAML(this_fpt)
        fpt.append(f)

    fpt_n = len(fpt)
    print(f'{fpt_n} fingerprints read')

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
