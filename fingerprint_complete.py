"""=================================================================================================
Generate a fingerprint from a XIOS file. XIOS structures are first broken into quasi-disjoint 
segments, then all possible subsets of vertices enumerated for each segment.

This code was modified from the old fingerprint random.py. Although this is the first version, it is
labeled as 2.2 to be consistent with the overall XIOS version

9 April 2026     Michael Gribskov
================================================================================================="""
version = 2.2

import sys
import os
import datetime
import argparse
from itertools import combinations
from collections import defaultdict
from topology import Topology
from xios import Gspan, MotifDB, Xios
from fingerprint import Fingerprint


def process_command_line():
    """---------------------------------------------------------------------------------------------
    usage: fingerprint_complete.py [-h] [-m MOTIFDB] [-r RNA] [-f FPT] [-s SUBGRAPHSIZE] [-c COVERAGE]
                                 [-l LIMIT] [-n] [-q]

    Calculate an XIOS fingerprint from a XIOS XML file by random sampling

    optional arguments:
      -h, --help                            show this help message and exit
      -m MOTIFDB, --motifdb MOTIFDB         RNA motif file
      -r RNA, --rna RNA                     RNA topology in XIOS XML format
      -f FPT, --fpt FPT                     Fingerprint output file (default=STDOUT)
      -s SUBGRAPHSIZE, --subgraphsize       SUBGRAPHSIZE
                                            Size subgraph to sample (default=7)
      -c COVERAGE, --coverage COVERAGE     Maximum segment size  (default=20)
      -n, --noparent                        Exclude parent graphs from fingerprint (default=False)
      -q, --quiet                           Minimal output on stdout (default=False)

    :return:
    ---------------------------------------------------------------------------------------------"""
    cl = argparse.ArgumentParser(
        description=f'v{version} Calculate XIOS fingerprint fromXIOS file by segmentation and complete sampling',
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
                    default=7)
    cl.add_argument('-c', '--coverage',
                    help='Maximum segment size (default=%(default)s)',
                    type=int,
                    default=20)
    cl.add_argument('-n', '--noparent',
                    help='Do not include parent graphs in fingerprint (default=%(default)s)',
                    action='store_true')
    cl.add_argument('-b', '--basesmin',
                    help='Minimum number of paired bases in a stem (%(default)s)',
                    default=2, type=int)
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

    # identify all edges between the vertices in vlist
    for r in range(len(vlist) - 1):
        row = vlist[r]
        for c in range(r + 1, len(vlist)):
            col = vlist[c]
            if adj[row][col] in 'ijo':
                struct.append([row, col, edge[adj[row][col]]])

    # print(struct)
    return Xios(list=struct)


def connected(adj, vlist):
    """---------------------------------------------------------------------------------------------
    True if vertices are connected

    :param adj: Xios.adjacency      adjacency matrix
    :param vlist: list              list of vertices to test
    :return: bool                   True if connected
    ---------------------------------------------------------------------------------------------"""
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
    """---------------------------------------------------------------------------------------------
    For each stem in the structure, construct a list of vertex sets, vset_list. vset_list[v] is the
    set of vertices
    that have non-serial relationships to v

    :param xios: XIOS       RNA structure
    :return: list           list of vsets
    ---------------------------------------------------------------------------------------------"""
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
    """---------------------------------------------------------------------------------------------
    extend a vertex set,
    for a source vertex v, the vertex set is vset_list[v]. the extended set is the union of 
    vset_list[v] with all vertex sets vset_list[u] that occur in vset_list[v] and u > v vertices in
    vset_list[v] with u < v are simply added to the merged set.
    
    Rationale:
    The stems (vertices) in the adjacency list are sorted by the coordinate of their leftmost 
    basepair. This means a vertex u with u<v is either unconnected, contains v (vju), or
    pseudoknotted (vou). If u is merged into the segment with union, it will include all the stems
    contained within the bigger outside stem (u). In this case the only segments found will be the
    completely disjoint subgraphs. A vertex w with v<w is either contained in v(viw), or
    pseudoknotted. W can be added to the segment using union.

    :param v: int             index of a source vertex set in vset_list
    :param vset_list: list    list of vsets for each vertex (defined by get_vset_list)
    :return: set              extended set of vertices
    ---------------------------------------------------------------------------------------------"""
    extended = vset_list[v]
    for u in list(vset_list[v]):
        x=1
        # extended.add(u)
        if u > v:
            extended = extended.union(vset_list[u])

    return extended


def vset_merge(candidate, limit):
    """---------------------------------------------------------------------------------------------
    vset_merge() identifies vertex sets that are quasi-disjoint. That is, they are densely connected
    inside the vset and loosely connected (quasi-disjoint) between vsets. These vsets are candidates
    for subgraphs that are small enough to enumerate all the k-vertex subgraphs, but that preserve
    most of the important information from the original graph.

    Vset_merge() compares candidate segments and identifies those that can be combined (union)
    without exceeding the size limit. Greedy approach: starting from the smallest set, add to all
    larger sets that have an intersection, subject to the size limit. When candidate vsets can no
    longer be merged (because they have no intersection with others, or because they have reached
    the size limit, they are moved to the final (approved) list.

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

        after checking all vertex sets in current, if no merges are found, add smallest to the final
        selected set

    :param candidate:list of set    vset_list, vertex subsets from segment_vset()
    :param limit: int               maximum number of elements in segment sets
    :return: list                   vset_list final list of vsets to segment structure
    ---------------------------------------------------------------------------------------------"""
    final = []
    current = sorted(candidate, key=lambda x: len(x), reverse=True)
    while current:
        small = current.pop()
        merged = False
        for v, s in enumerate(current):
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
    """---------------------------------------------------------------------------------------------
    Divide a structure into segments based on the sets of vertices that are iox connected to each
    vertex (vsets). This method detects all strictly disjoint sets of vertices, as well as 
    quasi-disjoint sets that have very few connecting edges to others.
    Method:
        make a list of all available vsets, process until this list is empty 
            (all vertices occur in a segment or the set is disjoint)
        process vertex lists in decreasing size order
            pop largest vertex set from the list (branches)
            if largest > limit, can't be a segment, get the next largest vertex set from the list
            if largest has no intersection with the available vertices, these vertex sets are
                already covered, next

            extend this vertex set: extended = vset_extend(largest_vset, branches)
            if extended > size limit:
                not a good segment, discard
            else:
                add to segment list
                remove segment from available vertices

        merge segments in segment list where possible: segment = vset_merge(segment, limit)

    :param xios: XIOS     RNA structure
    :param limit: int     maximum vertices, should be a size that can be exhaustively sampled
    :return: list         list of vertex sets, the selected quasi-disjoint vertex subsets
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

        # extended = vset_extend(bigv, branches)
        # if len(extended) > limit:
        #     continue

        segment.append(biggest)
        # remove from available
        available.difference_update(biggest)
        # print(available)

    segment = vset_merge(segment, limit)
    return segment


# ##################################################################################################
# Main
# ##################################################################################################
start = datetime.datetime.now()
runstart = start.strftime('%Y-%m-%d %H:%M:%S')
opt = process_command_line()

verbose = not opt.quiet
mode = None
if opt.noparent:
    mode = 'simple'
else:
    mode = 'extended'

if opt.quiet:
    print(f'fingerprint_complete: {runstart};motifdb:{opt.motifdb.name};{opt.rna.name};fpt:{opt.fpt}',
          end='\n')
else:
    print(f'fingerprint_complete v{version} - Calculate XIOS fingerprint from RNA topology {runstart}')
    print('\tRNA structure: {}'.format(opt.rna.name))
    print('\tMotif database: {}'.format(opt.motifdb.name))
    print('\tFingerprint: {}'.format(opt.fpt))
    print(f'\tMinimum stem size: {opt.basesmin}')
    print('\tSubgraph size: {}'.format(opt.subgraphsize))
    print(f'\tMaximum segment size: {opt.coverage}')
    if opt.noparent:
        print('\tMode: simple, parent motifs are omitted')
    else:
        print('\tMode: extended, parent motifs are included')

# read in the RNA structure, and remove stems with less than minimum number of basepairs
rna = Topology(xml=opt.rna)
xios_filter(rna, opt.basesmin)

fingerprint = Fingerprint()
fingerprint.information['Date'] = runstart
fingerprint.information['Run details'] = {'version': f'fingerprint_complete v{version}'}
fingerprint.information['Run details']['Maximum segment size'] = opt.coverage
fingerprint.information['Run details']['Mode'] = mode
fingerprint.information['Run details']['Minimum number of stems'] = opt.basesmin
fingerprint.information['Run details']['Maximum segment size'] = opt.coverage
fingerprint.information['Run details']['Subgraph size'] = opt.subgraphsize
# motifdb and xios file are opened by argparse, so you must use the name attribute
fingerprint.information['RNA structure'] = opt.rna.name
fingerprint.information['Fingerprint file'] = opt.fpt
# below, use .name because these files are opened by arg_parse
if not opt.noparent:
    # only use motifDB if parents are required
    fingerprint.information['Motif database'] = {'name': opt.motifdb.name}
fingerprint.information['Result'] = {'Number of stems after filtering': len(rna.stem_list)}

if verbose: print(f'\n\tSegmenting structure with {len(rna.stem_list)} stems')

selected = defaultdict(int)
adj = rna.adjacency
count = 0
# segment structure into tractable pieces
segment = segment_vset(rna, opt.coverage)
if verbose:
    plural = ''
    if len(segment) > 1: plural = 's'
    print(f'\t{len(segment)} segment{plural} found')
    for s in segment:
        print(f'\t{s}')

if verbose: print(f'\n\tIdentifying {opt.subgraphsize} stem motifs')
fingerprint.information['Result']['Segments'] = len(segment)

# for each segment > basemin (minimum number of stems) exhaustive enumerate combinations and convert
# to DFS code (motifs)
dfsset = defaultdict(int)
adjacency = rna.adjacency
nsample = 0
n_disjoint = 0
for subgraph in segment:
    subgraphsize = min(len(subgraph) - 1, opt.subgraphsize)
    subgraphsize = max(3, subgraphsize)
    if verbose: print(f'({subgraphsize})\tsubgraph:{sorted(subgraph)}')

    ssample = 0
    for vlist in combinations(list(subgraph), subgraphsize):
        # sample all combinations of size subgraphsize
        nsample += 1
        ssample += 1

        if not connected(adjacency, vlist):
            n_disjoint += 1
            # print(f'\t{n_disjoint} disconnected: {vlist}')
            continue

        if verbose and not nsample % 10000:
            print(f'\t\t{nsample}:{ssample}:{n_disjoint} {vlist}\t{subgraph}')

        xios = expand_vertices(vlist, adj)
        gspan = Gspan(graph=xios)
        dfsset[gspan.minDFS().human_encode()] += 1

    if verbose: print(f'\t\t{nsample}:{ssample}:{n_disjoint} motifs: {len(dfsset)}')

for dfs in sorted(dfsset):
    # finally store the fingerprints and counts
    fingerprint.add(dfs, dfsset[dfs])

# Summary and notes for fingerprint file
simple_n = fingerprint.n
fingerprint.information['Result']['Simple fingerprints'] = simple_n
fingerprint.information['Result']['Total subgraphs'] = nsample
fingerprint.information['Result']['Disjoint subgraphs'] = n_disjoint

if verbose:
    print(f'\n\tOverall sampling stats: total:{nsample}  disjoint:{n_disjoint}')
    print(f'\tSimple fingerprint motifs: {simple_n}')
else:
    if opt.noparent:
        # for quiet mode only print if no extended motifs
        print(f'\nsegment_size:{opt.coverage}|subgraph:{opt.subgraphsize}|simple:{simple_n}', end='')

# Generate parent motifs from motif database if requested

if opt.noparent:
    pass
else:
    # add the parents
    motif = MotifDB.unpickle(opt.motifdb)
    fingerprint.information['Motif database']['checksum'] = motif.information['checksum']
    fingerprint.information['Motif database']['description'] = motif.information['name']
    extended_n = fingerprint.add_parents(motif)
    fingerprint.information['Result']['Extended fingerprints'] = extended_n

    if opt.quiet:
        print(f'\nsegment_size:{opt.coverage}|subgraph:{opt.subgraphsize}|', end='')
        print(f'sample:{fingerprint.count}|simple:{simple_n}|extended:{fingerprint.n}')
    else:
        print(f'\tExtended fingerprint motifs:{fingerprint.n}')

if verbose: print(f'\twriting to {opt.fpt}')

fingerprint.writeYAML(opt.fpt)
stop = datetime.datetime.now()
runend = stop.strftime('%Y-%m-%d %H:%M:%S')
runtime = stop - start
fingerprint.information['Elapsed time'] = str(runtime)

if verbose:
    print('\nCompleted: {}'.format(runend))
    print(f'Elapsed time: {runtime}')

exit(0)
