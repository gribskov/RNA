"""=================================================================================================
sample a random fingerprint from a structure in XIOS XML format

10 July 2019     Michael Gribskov
================================================================================================="""
import sys
import os
import datetime
import argparse
import ast
from collections import defaultdict
from topology import Topology
from xios import Gspan, MotifDB, Xios
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
    default_subgraphsize = 7
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
    cl.add_argument('-b', '--basesmin',
                             help='Minimum number of paired bases in a stem (%(default)s)',
                             default=3, type=int)
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


def xios_filter(rna, minbp):
    """---------------------------------------------------------------------------------------------
    remove stems with less than minbp paired bases. lvienna shows base-pairs as ( so count the
    number of ( to get the number of basepairs

    :param rna: Xios        structrure to filter
    :param minbp: int       minimum number of base-pairs in stem
    :return: int            number of stems
    ---------------------------------------------------------------------------------------------"""
    for stem in list(rna.stem_list):
        if stem.lvienna.count('(') < minbp:
            rna.stem_list.remove(stem)

    return len(rna.stem_list)

def expand_vertices(vlist, adj):
    """---------------------------------------------------------------------------------------------
    given a set of vertices and adjacency matrix, select all edges and return the xios

    :param vlist: list          integers indicating the vertices in the structure graph
    :param adj: list of list    adjacency matrix
    :return: Xios
    ---------------------------------------------------------------------------------------------"""
    edge = {'i': 0, 'j': 1, 'o': 2, 's': 3, 'x': 4}
    struct = []

    # identify all the edges between the vertices in vlist
    for r in range(len(vlist) - 1):
        row = vlist[r]
        for c in range(r + 1, len(vlist)):
            col = vlist[c]
            if adj[row][col] in 'ijo':
                struct.append([row, col, edge[adj[row][col]]])

    # print(struct)
    return Xios(list=struct)

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
    print(f'\tMinimum stem size: {opt.basesmin}')

# read in the RNA structure
rna = Topology(xml=opt.rna)
xios_filter(rna, opt.basesmin)
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
fingerprint.information['Number of stems'] = len(rna.stem_list)
fingerprint.information['Minimum number of stems'] = opt.basesmin

if len(rna.stem_list) < opt.subgraphsize:
    opt.subgraphsize = max(len(rna.stem_list)-2,1)
    print(f'Reset subgraph size to {opt.subgraphsize}')

from itertools import combinations
from  math import comb
stemlist = rna.stem_list
expect =  comb(len(rna.stem_list),opt.subgraphsize)
print(f'expect {expect} sets from {len(rna.stem_list)} stems')
excluded = 0
included = 0

print(f'\n\tGenerating connected vertex sets')
selected = defaultdict(int)
adj = rna.adjacency
count = 0

# different vertex selections may produce the same dfs so collect them here
print(f'\nCollecting DFS from vertex sets ({len(selected)})')
dfsset = defaultdict(int)
stemlist = rna.stem_list
expect =  comb(len(rna.stem_list),opt.subgraphsize)
print(f'expect {expect} sets from {len(rna.stem_list)} stems')
count = 0
for vlist in combinations(range(len(rna.stem_list)), opt.subgraphsize):
    count += 1
    # vset = ast.literal_eval(vstr)
    begin = stemlist[vlist[0]].lbegin
    end = stemlist[vlist[0]].rend
    disjoint = False
    for v in vlist:
        stem = stemlist[v]
        if stem.lbegin < end:
            end = max(end, stem.rend)
        else:
            disjoint = True
            break

    if not disjoint:
        print(vlist)
        xios = expand_vertices(vlist, adj)
        gspan = Gspan(graph=xios)
        dfsset[gspan.minDFS().human_encode()] += 1
        dfs = gspan.minDFS().human_encode()
        dfsset[dfs] += 1
        if not count % 1000:
            print(f'{count}\t{dfs}\t{dfsset[dfs]}')
    else:
        print(f'vlist ({vlist}) is disjoint')

print(f'Uniquifying DFS codes')
for dfs in dfsset:
    # finally store the fingerprints and counts
    fingerprint.add(dfs, dfsset[dfs])

if opt.quiet:
    print(f'\ncov:{opt.coverage}|subgraph:{opt.subgraphsize}|limit:{opt.limit}', end='')
else:
    fingerprint.information['Coverage'] = opt.coverage
    fingerprint.information['Subgraph size'] = opt.subgraphsize
    fingerprint.information['Sampling limit'] = opt.limit
    fingerprint.information['Sample size'] = fingerprint.count
    print(f'\tSample size: {fingerprint.count}')

# to include parent, you must read a motif database.  this is only done after all the motifs have
# been added to the fingerprint

simple_n = fingerprint.n
if opt.noparent:
    pass
else:
    # add the parents
    motif = MotifDB.unpickle(opt.motifdb)
    extended_n = fingerprint.add_parents(motif)
    fingerprint.information['Motif database checksum'] = motif.information['checksum']
    fingerprint.information['Motif database description'] = motif.information['name']
    fingerprint.information['Motif database checksum'] = motif.information['checksum']

    if opt.quiet:
        print(f'|sample:{fingerprint.count}|xpt:{fingerprint.n}', end='')
    else:
        print(f'\tExtended fingerprint motifs:{fingerprint.n}')

if opt.quiet:
    print(f'|fpt:{simple_n}')
else:
    # write simple fpt information regardless of whether parents are added
    print(f'\tSimple fingerprint motifs: {simple_n}')
    print(f'\twriting to {opt.fpt}')

fingerprint.writeYAML(opt.fpt)
daytime = datetime.datetime.now()
runend = daytime.strftime('%Y-%m-%d %H:%M:%S')
sys.stdout.write('Completed: {}'.format(runend))

exit(0)
