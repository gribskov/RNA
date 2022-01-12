"""=================================================================================================
sample a random fingerprint from a structure in XIOS XML format

10 July 2019     Michael Gribskov
================================================================================================="""
import sys
import datetime
import argparse
from topology import Topology
from xios import Xios, Gspan, MotifDB
from fingerprint import Fingerprint


def process_command_line():
    """---------------------------------------------------------------------------------------------
    usage: fingerprint_random.py [-h] [--motifdb_file MOTIFDB_FILE] [--rna_file RNA_FILE] [--fpt_file FPT_FILE]
                                 [--subgraphsize SUBGRAPHSIZE] [--coverage COVERAGE] [--noparent]

    Calculate an XIOS fingerprint from a XIOS XML file by random sampling

    optional arguments:
      -h, --help                   show this help message and exit
      --motifdb_file MOTIFDB_FILE  RNA motif file
      --rna_file RNA_FILE          RNA topology in XIOS XML format
      --fpt_file FPT_FILE          Fingerprint output file (default=STDOUT)
      --subgraphsize SUBGRAPHSIZE  Size subgraph to sample (default=6)
      --coverage COVERAGE          Minimum coverage for sampled graphs (default=3)
      --noparent                   Include parent graphs in fingerprint (default=False)

    :return:
    ---------------------------------------------------------------------------------------------"""
    default_subgraphsize = 6
    default_coverage = 3

    cl = argparse.ArgumentParser(
        # argument_default=argparse.SUPPRESS,
        # formatter_class=lambda prog:argparse.ArgumentDefaultsHelpFormatter(prog,width=200,max_help_position=40),
        description='Calculate an XIOS fingerprint from a XIOS XML file by random sampling',
        formatter_class=lambda prog: argparse.HelpFormatter(prog,width=120,max_help_position=40)
    )
    cl.add_argument('--motifdb',
                    help='RNA motif file',
                    type=argparse.FileType('r'))
    cl.add_argument('--rna',
                    help='RNA topology in XIOS XML format',
                    type=argparse.FileType('r'))
    cl.add_argument('--fpt',
                    help='Fingerprint output file (default=STDOUT)',
                    type=argparse.FileType('w'),
                    default=sys.stdout)
    cl.add_argument('--subgraphsize',
                    help='Size subgraph to sample (default=%(default)s)',
                    type=int,
                    default = default_subgraphsize)
    cl.add_argument('--coverage',
                    help='Minimum coverage for sampled graphs (default=%(default)s)',
                    type=int,
                    default=default_coverage)
    cl.add_argument('--noparent',
                    help='Include parent graphs in fingerprint (default=%(default)s)',
                    action='store_true')

    return cl.parse_args()


# ##################################################################################################
# Main
# ##################################################################################################
daytime = datetime.datetime.now()
runstart = daytime.strftime('%Y-%m-%d %H:%M:%S')
opt = process_command_line()

# for testing
# opt.motifdb_file = 'data/5stem.list.txt'
# opt.rna_file = 'data/rnasep_a1.Buchnera_APS.xios'
# opt.fpt_file = 'data/test2.fpt'
print('fingerprint_random - Sample XIOS fingerprint from RNA topology {}'.format(runstart))
print('\tRNA structure: {}'.format(opt.rna.name))
# print('\tMotif database: {}'.format(opt.motifdb.name))
print('\tFingerprint: {}'.format(opt.fpt.name))
print('\tSubgraph size: {}'.format(opt.subgraphsize))
print('\tCoverage (minimum): {}'.format(opt.coverage))
print('\tOmit parents: {}'.format(opt.noparent))

# read in the motif database and RNA structure
# motif = MotifDB(json=opt.motifdb)
rna = Topology(xml=opt.rna)

# this is an unweighted sampling strategy.  Others were tried, sampling:
#   inversely proportional to number of times previously sampled
#       scaled by 1/n and 1/rank
#   proportional to number of neighbors
#   inversely proportional to number of neighbors

fingerprint = Fingerprint()
fingerprint.information['Date'] = runstart
fingerprint.information['File'] = opt.fpt.name
# fingerprint.information['Motif database'] = opt.motifdb.name
# fingerprint.information['Motif database checksum'] = motif.information['checksum']
# fingerprint.information['Motif database description'] = motif.information['name']
fingerprint.information['RNA structure'] = opt.rna.name

subgraph = opt.subgraphsize
count_threshold = opt.coverage

# sample the first subgraph to have a motif with minimum occurance
xios = rna.sample_xios(subgraph)
gspan = Gspan(graph=xios)
dfs = gspan.minDFS().human_encode()
print( dfs)
fingerprint.add(dfs)
minmotif = fingerprint.minkey()

while True:
    # sample until the lowest count motif is above the opt.coverage count_threshold.  You only have
    # to recheck the minimum count when your current minimum graph passes the threshold (finding
    # minimum is expensive)
    xios = rna.sample_xios(subgraph)
    if len(xios) < subgraph:
        # in case the sampled graph is empty
        continue

    gspan = Gspan(graph=xios)
    dfs = gspan.minDFS().human_encode()
    if not fingerprint.total % 10000:
        print(fingerprint.total, dfs)
        fingerprint.writeYAML(sys.stderr)
    fingerprint.add(dfs)

    if dfs == minmotif:
        minmotif = fingerprint.minkey()
        mincount = fingerprint.mincount()
        print('{}\t{}\t{}'.format(fingerprint.total, fingerprint.nmotif, fingerprint.mincount()))
        if mincount >= count_threshold:
            break

print('Simple fingerprint: {}\t{}\t{}'.format(fingerprint.total, fingerprint.nmotif,
                                              fingerprint.mincount()))
if not opt.noparent:
    fingerprint.add_parents(motif)
    print('Extended fingerprint: {}\t{}\t{}'.format(fingerprint.total, fingerprint.nmotif,
                                                    fingerprint.mincount()))

# print(fingerprint.toYAML())
fingerprint.writeYAML(opt.fpt)

daytime = datetime.datetime.now()
runend = daytime.strftime('%Y-%m-%d %H:%M:%S')
sys.stdout.write('Completed: {}'.format(runend))

exit(0)
