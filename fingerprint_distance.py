"""=================================================================================================
Calculate the distance between fingerprints.
Assume fingerprint are all in the same directory
Calculate the jaccard distance between all pairs

Michael Gribskov     04 February 2022
================================================================================================="""
import glob
from os.path import basename
from fingerprint import Fingerprint, FingerprintSet
def process_command_line():
    """---------------------------------------------------------------------------------------------

    ---------------------------------------------------------------------------------------------"""
    cl = argparse.ArgumentParser(
        description='Calculate distances between fingerprints',
        formatter_class=lambda prog: argparse.HelpFormatter(prog, width=120, max_help_position=40)
    )
    cl.add_argument('-m', '--dir',
                    help='Directory with fingerprint files (default=%(default)s)',
                    default='../data/fpt')
    cl.add_argument('-r', '--suffix ',
                    help='Suffix for fingertprint files  (default=%(default)s',
                    tdefault='.fpt')
    cl.add_argument('-f', '--fpt',
                    help='Fingerprint output file (default=STDOUT, auto=calculate from xios)',
                    default='auto')
    cl.add_argument('-o', '--out',
                    help='Output file name (default=%(default)s)',
                    default='distance.out')


    args = cl.parse_args()

    if str(args.fpt) == 'auto':
        args.fpt = fpt_from_xios(args)

    return args
# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    daytime = datetime.datetime.now()

    opt = process_command_line()
    fpt_dir = opt.fpt'
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
    # for i in range(fpt_n):
    #     fpt1 = fpt[i].information['File']
    #     for j in range(i + 1, fpt_n):
    #         fpt2 = fpt[j].information['File']
    #         j = fpt.jaccard_sim([i,j])
    #         bc = fpt.bray_curtis_dis([i,j])
    #         print(f'{fpt1}\t{fpt2}\t{j[2]:.3f}\t{bc:[2]:.3f}')
    maximum = {'jaccard':0, 'bray-curtis':0}
    minimum = {'jaccard': 1000000, 'bray-curtis': 1000000}
    
    jaccard = fpt.jaccard_sim([])
    bc = fpt.bray_curtis_dis([])

    out = open( outfile, 'w')
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
