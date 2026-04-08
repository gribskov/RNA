"""=================================================================================================
sample a random fingerprint from a structure in XIOS XML format

10 July 2019     Michael Gribskov
================================================================================================="""
import sys
import os
import datetime
import argparse
import ast
from itertools import combinations
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

    :param rna: Xios        structure to filter
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

def connected(adj, vlist):
    """-----------------------------------------------------------------------------------------------------------------
    True if vertices are connected

    :param adj: Xios.adjacency      adjacency matrix
    :param vlist: list              list of vertices to test
    :return: bool                   True if connected
    -----------------------------------------------------------------------------------------------------------------"""
    # if len(subset) <= 1:
    #     return True

    subset_set = set(vlist)
    visited = {vlist[0]}
    stack = [vlist[0]]

    while stack:
        u = stack.pop()
        # iterate only over vertices in subset
        for v in vlist:
            if v not in visited and adj[u][v] != 's':
                visited.add(v)
                stack.append(v)

            if len(visited) == len(vlist):
                return True

    return False



def get_vset_list(xios):
    """-----------------------------------------------------------------------------------------------------------------
    For each stem in the structure, construct a list of vertex sets, vset_list. vset_list[v] is the set of vertices
    that have non-serial relationships to v

    :param xios: XIOS       RNA structure
    :return: list           list of vsets
    -----------------------------------------------------------------------------------------------------------------"""
    adj = xios.adjacency
    vset_list = []
    for i in range(len(adj)):
        vset_list.append(set())
        for j in range(len(adj[i])):
            # if adj[i][j] in 's-':
            if adj[i][j] == 's':
                # do not include in vset
                continue
            vset_list[-1].add(j)

    return vset_list


def vset_extend(v, vset_list):
    """-----------------------------------------------------------------------------------------------------------------
    extend a vertex set,
    for a source vertex v, the vertex set is vset_list[v]
    the extended set is the union of vset_list[v] with all vertex sets vset_list[u] that occur in vset_list[v] and u > v
    vertices in vset_list[v] with u < v are simple added to the merged set

    :param v: int           index of a source vertex set in vset_list
    :param vset_list: list   list of sets of vertices connected to each vertex (defined by stem_overlap)
    :return: set            extended set of vertices
    -----------------------------------------------------------------------------------------------------------------"""
    extended = vset_list[v]
    for thisv in list(vset_list[v]):
        if thisv > v:
            extended = extended.union(vset_list[v])
        else:
            extended.add(thisv)

    return extended


def vset_merge(candidate, limit):
    """-----------------------------------------------------------------------------------------------------------------
    segment_vset() identifies vertex sets that are quasi-disjoint. That is, they are densely connected inside the vset
    and loosely connected or disjoint between vsets. These vsets are candidates for subgraphs that are small enough to
    enumerate all the k-vertex subgraphs, but that preserve most of the important information from the original graph.

    Vset_merge() compares candidate vsets and identifies those that can be combined (union) without exceeding
    the size limit. Greedy approach: starting from the smallest set, add to all larger sets that have an intersection,
    subject to the size limit. When candidate vsets can no longer be merged (because they have no intersection with
    others, or because they have reached the size limit, they are moved to the final (approved) list.

    sort candidate vset list (candidate) by size, push on current stack
    add all vsets to the list of current vsets
    working from the smallest vertex set to the largest:
        pop the smallest vertex set from the current stack
        for each vertex set remaining in current:
            calculate intersection with smallest
            if not intersection, continue

            calculate union with smallest
            if union > size_limit, continue

            otherwise, smallest can extend the current vertex set
                replace current with union
                remove smallest from available vertices

        after checking all vertex sets in current, if no merges are found, add smallest to the final selected set

    :param candidate:list of set    vset_list, vertex subsets from segment_vset()
    :param limit: int               maximum number of elements in segment sets
    :return: list                   vset_list final list of vsets to segment structure
    -----------------------------------------------------------------------------------------------------------------"""
    final = []
    current = sorted(candidate, key=lambda x: len(x), reverse=True)
    while current:
        small = current.pop()
        merged = False
        for v,s in enumerate(current):
            if small.intersection(s):
                # overlap, try to combine
                testmerge = small.union(s)
                if len(testmerge) <= limit:
                    current[v] = testmerge
                    merged = True

        if not merged:
            # if small cannot be merged it must be saved
            final.append(small)

    return final


def segment_vset(xios, limit):
    """-----------------------------------------------------------------------------------------------------------------
    Divide a structure into segments based on the sets of vertices that are iox connected to each vertex (i.e.,
    branches). This method detects all strictly disjoint sets of vertices so no other check is needed.
    Method:
        make a list of all available vertices, process until this list is empty (all vertices occur in a segment)
        process vertex lists in decreasing size order
            pop largest vertex set from the list (branches)
            if largest > limit, can't be a segment, get the next largest vertex set from the list
            if largest has no intersection with the available vertices, these vertex sets are already covered, next

            extend this vertex set: extended = vset_extend(largest_vset, branches)
            if extended > size limit:
                not a good segment, discard
            else:
                add to segment list
                remove segment from available vertices

        merge segments in segment list where possible: segment = vset_merge(segment, limit)

    :param xios: XIOS     RNA structure
    :param limit: int     maximum number of vertices, this should be a size that can be exhaustively sampled
    :return: list         list of vertex sets, the final set of selected quasi-disjoint vertex subsets
    -----------------------------------------------------------------------------------------------------------------"""
    segment = []
    available = set([i for i in range(len(xios.adjacency))])
    branches = get_vset_list(xios)

    # sort indices to the branches list (list of vertex sets) by size of set
    order = [i for i in range(len(branches))]
    sorder = sorted(order, key=lambda x: len(branches[x]))

    while available:
        try:
            bigv = sorder.pop()
            biggest = branches[bigv]
        except IndexError:
            break
        if len(biggest) > limit:
            continue

        if not biggest.intersection(available):
            # there are no unique vertices here
            continue

        extended = vset_extend(bigv, branches)
        if len(extended) > limit:
            continue

        segment.append(extended)
        # remove from available
        available.difference_update(extended)
        # print(available)

    segment = vset_merge(segment, limit)
    return segment

# ##################################################################################################
# Main
# ##################################################################################################
daytime = datetime.datetime.now()
runstart = daytime.strftime('%Y-%m-%d %H:%M:%S')
opt = process_command_line()
segment_limit = 25

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
    print(f'\tSegment limit: {segment_limit}')

# read in the RNA structure
rna = Topology(xml=opt.rna)
xios_filter(rna, opt.basesmin)
# print(rna.format_edge_list())

fingerprint = Fingerprint()
fingerprint.information['Date'] = runstart
fingerprint.information['File'] = opt.fpt
# below, use .name because these files are opened by arg_parse
fingerprint.information['Motif database'] = opt.motifdb.name
fingerprint.information['RNA structure'] = opt.rna.name
fingerprint.information['Number of stems'] = len(rna.stem_list)
fingerprint.information['Minimum number of stems'] = opt.basesmin

print(f'\n\tGenerating connected vertex sets')
selected = defaultdict(int)
adj = rna.adjacency
count = 0

# segment structure into tractable pieces
segment = segment_vset(rna, 25)

# for each segment > basemin (minimum number of stems) exhaustive enumerate combinations and convert
# to DFS code (motifs)
nsample = 0
n_disjoint = 0
dfsset = defaultdict(int)
adjacency = rna.adjacency
for subgraph in segment:
    print(f'subgraph:{subgraph}')
    ssample = 0
    for vlist in combinations(list(subgraph),opt.subgraphsize):
        # sample all combination of size opt.limit
        nsample += 1
        ssample += 1

        if not connected(adjacency, vlist):
            n_disjoint += 1
            # print(f'\t{n_disjoint} disconnected: {vlist}')
            continue

        if not nsample % 10000:
            print(f'\t{nsample}:{ssample}:{n_disjoint} {vlist}\t{subgraph}')

        xios = expand_vertices(vlist, adj)
        gspan = Gspan(graph=xios)
        # print(gspan.graph)
        dfsset[gspan.minDFS().human_encode()] += 1

    print(f'\t{nsample}:{ssample}:{n_disjoint} motifs: {len(dfsset)}')

# modified the code to first generate samples of connected vertices (vertex sets) and only later
# generate the DFS codes. This saves a lot of time for small graphs, but not as much for large
# graphs where the number of combinations is much larger than the sample so most vertex sets are unique
# TODO one possibility would be to convert to DFS in batches
# for _ in range(opt.limit):
#     # sample until the lowest count motif is above the opt.coverage count_threshold.  You only have
#     # to recheck the minimum count when your current minimum graph passes the threshold (finding
#     # minimum is expensive)
#     xios = rna.sample_all(adj, opt.subgraphsize)
#     if not xios:
#         # should never reach here, revised sample_xios will always return a graph
#         sys.stderr.write(f'fingerprint_random - graph could not be sampled, possibly too small')
#         exit(2)
#
#     selected[str(xios)] += 1
#     # print(f'min:{min(selected.values())}')
#     minselect = min(selected.values())
#     count += 1
#     if not count % 10000:
#         print(f'\t\t{count} sets\t minimum count {minselect}')
#     if minselect >= opt.coverage:
#         # TODO when saving vertex lists,  this test will never me met
#         break
#
# # different vertex selections may produce the same dfs so collect them here
# print(f'\nCollecting DFS from vertex sets ({len(selected)})')
# dfsset = defaultdict(int)
# minsubgraph = 4
# for vstr in selected:
#     vset = ast.literal_eval(vstr)
#     if len(vset) < minsubgraph:
#         continue
#     xios = expand_vertices(vset, adj)
#     gspan = Gspan(graph=xios)
#     dfsset[gspan.minDFS().human_encode()] += selected[vstr]
#
# print(f'Uniquifying DFS codes')
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
