"""=================================================================================================
sample a random fingerprint from a structure in XIOS XML format

10 July 2019     Michael Gribskov
================================================================================================="""
import sys
import os
import datetime
import argparse
from topology import Topology
from xios import Xios, Gspan, MotifDB
from fingerprint import Fingerprint


def process_command_line():
    """---------------------------------------------------------------------------------------------
    usage: fingerprint_random.py [-h] [-m MOTIFDB] [-r RNA] [-f FPT] [-s SUBGRAPHSIZE] [-c COVERAGE] [-l LIMIT] [-n] [-q]

    Calculate an XIOS fingerprint from a XIOS XML file by random sampling

    optional arguments:
      -h, --help                            show this help message and exit
      -m MOTIFDB, --motifdb MOTIFDB         RNA motif file
      -r RNA, --rna RNA                     RNA topology in XIOS XML format
      -f FPT, --fpt FPT                     Fingerprint output file (default=STDOUT)
      -s SUBGRAPHSIZE, --subgraphsize SUBGRAPHSIZE
                                            Size subgraph to sample (default=6)
      -c COVERAGE, --coverage COVERAGE      Minimum coverage for sampled graphs (default=3)
      -l LIMIT, --limit LIMIT               Maximum number of random graphs to sample (default=10000)
      -n, --noparent                        Do not include parent graphs in fingerprint (default=False)
      -q, --quiet                           Minimal output on stdout (default=False)


    :return:
    ---------------------------------------------------------------------------------------------"""
    default_subgraphsize = 6
    default_coverage = 3
    default_sampling_limit = 10000

    cl = argparse.ArgumentParser(
        # argument_default=argparse.SUPPRESS,
        # formatter_class=lambda prog:argparse.ArgumentDefaultsHelpFormatter(prog,width=200,max_help_position=40),
        description='Calculate an XIOS fingerprint from a XIOS XML file by random sampling',
        formatter_class=lambda prog: argparse.HelpFormatter(prog, width=120, max_help_position=40)
        )
    cl.add_argument('-m', '--motifdb',
                    help='pickled RNA motif file (default=%(default)s)',
                    type=argparse.FileType('rb'),
                    default='../RNA/data/2to7stem.mdb.pkl')
    cl.add_argument('-r', '--rna',
                    help='RNA topology in XIOS XML format',
                    type=argparse.FileType('r'))
    cl.add_argument('-f', '--fpt',
                    help='Fingerprint output file (default=STDOUT, auto=calculate from xios)',
                    default='auto')
    cl.add_argument('-o', '--outputdir',
                    help='Directory for fingerprint output (default=%(default)s)',
                    default='./')
    cl.add_argument('-s', '--subgraphsize',
                    help='Size subgraph to sample (default=%(default)s)',
                    type=int,
                    default=default_subgraphsize)
    cl.add_argument('-c', '--coverage',
                    help='Minimum coverage for sampled graphs (default=%(default)s)',
                    type=int,
                    default=default_coverage)
    cl.add_argument('-l', '--limit',
                    help='Maximum number of random graphs to sample (default=%(default)s)',
                    type=int,
                    default=default_sampling_limit)
    cl.add_argument('-n', '--noparent',
                    help='Do not include parent graphs in fingerprint (default=%(default)s)',
                    action='store_true')
    cl.add_argument('-q', '--quiet',
                    help='Minimal output on stdout (default=%(default)s)',
                    action='store_true')

    args = cl.parse_args()

    print(f'auto={args.fpt}')
    if str(args.fpt) == 'auto':
        args.fpt = fpt_from_xios(args)

    return args


def fpt_from_xios(args):
    '''---------------------------------------------------------------------------------------------
    create a name for the output file by
        1) removing directory path
        2) changing the suffix .xios to .fpt
        
    :param xiosfile: string - input xios file
    :return: string - name for fingerprint file
    ---------------------------------------------------------------------------------------------'''
    fpt = os.path.basename(args.rna.name)
    l = len(fpt)
    if fpt.rindex('.xios') == l - 5:
        fpt = fpt[:-5]

    fpt += '.fpt'

    # if specified, add the output directory
    if args.outputdir:
        fpt = args.outputdir + fpt

    return fpt


# ##################################################################################################
# Main
# ##################################################################################################
daytime = datetime.datetime.now()
runstart = daytime.strftime('%Y-%m-%d %H:%M:%S')
opt = process_command_line()
print(f'motifdb:{opt.motifdb.name}')

if opt.quiet:
    print(f'fingerprint_random: {daytime} : {opt.rna.name} : {opt.fpt}', end=' ')
else:
    print('fingerprint_random - Sample XIOS fingerprint from RNA topology {}'.format(runstart))
    print('\tRNA structure: {}'.format(opt.rna.name))
    print('\tMotif database: {}'.format(opt.motifdb.name))
    print('\tFingerprint: {}'.format(opt.fpt.name))
    print('\tSubgraph size: {}'.format(opt.subgraphsize))
    print('\tCoverage (minimum): {}'.format(opt.coverage))
    print(f'\tMaximum sample: {opt.limit}')
    print('\tOmit parents: {}'.format(opt.noparent))

# read in the RNA structure
rna = Topology(xml=opt.rna)
# print(rna.format_edge_list())

# this is an unweighted sampling strategy.  Others were tried, sampling:
# inversely proportional to number of times previously sampled scaled by 1/n and 1/rank
# proportional to number of neighbors
# inversely proportional to number of neighbors

fingerprint = Fingerprint()
fingerprint.information['Date'] = runstart
fingerprint.information['File'] = opt.fpt
fingerprint.information['Motif database'] = opt.motifdb.name
fingerprint.information['RNA structure'] = opt.rna.name

# sample the first subgraph to have a motif with minimum occurence
xios = rna.sample_xios(opt.subgraphsize)
gspan = Gspan(graph=xios)
dfs = gspan.minDFS().human_encode()
fingerprint.add(dfs)
minmotif = fingerprint.minkey()

while True:
    # sample until the lowest count motif is above the opt.coverage count_threshold.  You only have
    # to recheck the minimum count when your current minimum graph passes the threshold (finding
    # minimum is expensive)
    xios = rna.sample_xios(opt.subgraphsize)

    gspan = Gspan(graph=xios)
    dfs = gspan.minDFS().human_encode()
    # print(dfs)
    if not opt.quiet and not fingerprint.count % 10000:
        print(fingerprint.count, dfs)
        # fingerprint.writeYAML(sys.stderr)
    fingerprint.add(dfs)

    if dfs == minmotif:
        # if the new dfs is the one with the lowest count, update the lowest count
        minmotif = fingerprint.minkey()
        mincount = fingerprint.mincount()
        # print('{}\t{}\t{}'.format(fingerprint.count, fingerprint.n, fingerprint.mincount()))
        if mincount >= opt.coverage or fingerprint.count > opt.limit:
            break

n_initial = fingerprint.n
if opt.noparent:
    if opt.quiet:
        print(f'fpt: {fingerprint.n} : {fingerprint.count}')
    else:
        print('Simple fingerprint: {}\t{}\t{}'.format(fingerprint.count, fingerprint.n, fingerprint.mincount()))

# to include parent, you must read a motif database
if not opt.noparent:
    # read in the motif database and RNA structure
    motif = MotifDB.unpickle(opt.motifdb)
    fingerprint.information['Motif database checksum'] = motif.information['checksum']
    fingerprint.information['Motif database description'] = motif.information['name']
    fingerprint.add_parents(motif)
    if opt.quiet:
        print(f'xpt: {fingerprint.n} : {fingerprint.count}')
    else:
        print('\nExtended fingerprint: {}\t{}\t{}'.format(fingerprint.count, fingerprint.n, fingerprint.mincount()))
        print(f'\t{n_initial} simple fingerprints extended to {fingerprint.n}')

# print(fingerprint.toYAML())
if not opt.quiet:
    print(f"\twriting to {motif.information['file']}")
print()
fingerprint.writeYAML(opt.fpt)

daytime = datetime.datetime.now()
runend = daytime.strftime('%Y-%m-%d %H:%M:%S')
sys.stdout.write('\nCompleted: {}'.format(runend))

exit(0)
