"""=================================================================================================
sample a random fingerprint from a structure in XIOS XML format

10 July 2019     Michael Gribskov
================================================================================================="""
import sys
import os
import datetime
import argparse
from topology import Topology
from xios import Gspan, MotifDB
from fingerprint import Fingerprint


def process_command_line():
    """---------------------------------------------------------------------------------------------
    usage: fingerprint_random.py [-h] [-m MOTIFDB] [-r RNA] [-f FPT] [-s SUBGRAPHSIZE] [-c COVERAGE]
                                 [-l LIMIT] [-n] [-q]

    Calculate an XIOS fingerprint from a XIOS XML file by random sampling

    optional arguments:
      -h, --help                            show this help message and exit
      -m MOTIFDB, --motifdb MOTIFDB         RNA motif file
      -r RNA, --rna RNA                     RNA topology in XIOS XML format
      -f FPT, --fpt FPT                     Fingerprint output file (default=STDOUT)
      -s SUBGRAPHSIZE, --subgraphsize       SUBGRAPHSIZE
                                            Size subgraph to sample (default=6)
      -c COVERAGE, --coverage COVERAGE      Minimum coverage for sampled graphs (default=3)
      -l LIMIT, --limit LIMIT               Maximum random graphs to sample (default=10000)
      -n, --noparent                        Exclude parent graphs from fingerprint (default=False)
      -q, --quiet                           Minimal output on stdout (default=False)


    :return:
    ---------------------------------------------------------------------------------------------"""
    default_subgraphsize = 6
    default_coverage = 3
    default_sampling_limit = 10000

    cl = argparse.ArgumentParser(
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

    if str(args.fpt) == 'auto':
        args.fpt = fpt_from_xios(args)

    return args


def fpt_from_xios(args):
    """---------------------------------------------------------------------------------------------
    create a name for the output file by
        1) removing directory path
        2) changing the suffix .xios to .fpt
        
    :param args: argparse namespace
    :return: string - name for fingerprint file
    ---------------------------------------------------------------------------------------------"""
    fpt = os.path.basename(args.rna.name)
    fptlen = len(fpt)
    if fpt.rindex('.xios') == fptlen - 5:
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

if opt.quiet:
    print(f'fingerprint_random: {daytime};motifdb:{opt.motifdb.name};{opt.rna.name};fpt:{opt.fpt}',
          end='\n')
else:
    print('fingerprint_random - Sample XIOS fingerprint from RNA topology {}'.format(runstart))
    print('\tRNA structure: {}'.format(opt.rna.name))
    print('\tMotif database: {}'.format(opt.motifdb.name))
    print('\tFingerprint: {}'.format(opt.fpt))
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
# below, use .name because these files are opened by arg_parse
fingerprint.information['Motif database'] = opt.motifdb.name
fingerprint.information['RNA structure'] = opt.rna.name

minmotif = ''
mincount = 0
while True:
    # sample until the lowest count motif is above the opt.coverage count_threshold.  You only have
    # to recheck the minimum count when your current minimum graph passes the threshold (finding
    # minimum is expensive)
    xios = rna.sample_xios(opt.subgraphsize)
    if not xios:
        print(f'graph could not be sampled, possibly too small')
        exit(2)
    gspan = Gspan(graph=xios)
    dfs = gspan.minDFS().human_encode()
    fingerprint.add(dfs)
    # print(dfs)

    if not opt.quiet and not fingerprint.count % 10000:
        # screen trace
        print(fingerprint.count, dfs)
        # fingerprint.writeYAML(sys.stderr)

    if (dfs == minmotif) or (not mincount):
        # if the new dfs is the one with the lowest count, update the lowest count, otherwise you
        # don't need to check
        minmotif = fingerprint.minkey()
        mincount = fingerprint.mincount()

    if mincount >= opt.coverage or fingerprint.count > opt.limit:
        # this is the successful exit point for the loop
        break

# to include parent, you must read a motif database.  this is only done after all the motifs have
# been added to the fingerprint
if opt.noparent:
    if opt.quiet:
        print(f'fpt: {fingerprint.n} : {fingerprint.count}')
    else:
        print('\tSimple fingerprint: {}\t{}\t{}'.format(fingerprint.count, fingerprint.n,
                                                      fingerprint.mincount()))
else:
    # add the parents
    motif = MotifDB.unpickle(opt.motifdb)
    simple_n = fingerprint.n
    extended_n = fingerprint.add_parents(motif)

    fingerprint.information['Motif database checksum'] = motif.information['checksum']
    fingerprint.information['Motif database description'] = motif.information['name']
    if opt.quiet:
        print(f'xpt: {fingerprint.n} : {fingerprint.count}')
    else:
        print(f'\nExtended fingerprint: {fingerprint.count}', end='\t')
        print(f'{fingerprint.n}', end='\t'),
        print(f'{fingerprint.mincount()}')
        print(f'{simple_n} simple fingerprints extended to {extended_n}')

if not opt.quiet:
    print(f'\twriting to {opt.fpt}')

fingerprint.writeYAML(opt.fpt)
daytime = datetime.datetime.now()
runend = daytime.strftime('%Y-%m-%d %H:%M:%S')
sys.stdout.write('Completed: {}'.format(runend))

exit(0)
