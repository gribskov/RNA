"""=====================================================================================================================
segment.py

Try some ideas for breaking large structures into parts

2026-04-02 gribskov
====================================================================================================================="""
import sys
import glob
from topology import Topology


def segment_strict(xios):
    """-----------------------------------------------------------------------------------------------------------------
    Check for completely disjoint substructures

    :param xios: Xios     RNA structure
    :return:
    -----------------------------------------------------------------------------------------------------------------"""
    stemlist = xios.stem_list
    segment = []
    begin = 0
    begin_pos = stemlist[begin].lbegin
    end = 0
    end_pos = stemlist[end].rend
    disjoint = False
    for i in range(1, len(stemlist)):
        next_begin = stemlist[i].lbegin
        next_end = stemlist[i].rend
        if next_begin < end_pos:
            # overlaps previous, update end
            end = i
            end_pos = max(end_pos, next_end)

        else:
            # does not overlap
            disjoint = True
            segment.append([begin, end])
            begin = end = i
            begin_pos = next_begin
            end_pos = next_end

    segment.append([begin, end])
    return segment


def segment_mountain(xios):
    """-----------------------------------------------------------------------------------------------------------------
    use the stem positions from lbegin to rend to construct a mountain plot

    :param xios: Xios     RNA structure
    :return:
    -----------------------------------------------------------------------------------------------------------------"""
    stemlist = xios.stem_list
    segment = []
    position = []
    n = 0
    for s in stemlist:
        position.append([s.lbegin, s.rend, n, 1])
        position.append([s.rend, s.lbegin, n, -1])
        n += 1

    height = 0
    n = 0
    for pos in sorted(position, key=lambda x: x[0]):
        height += pos[3]
        print(f'{n}\t{height}\t{pos}')
        n += 1

    return


def segment_line(xios):
    """-----------------------------------------------------------------------------------------------------------------
    select stems that cross a coordinate position

    :param xios: XIOS     RNA structure
    :return:
    -----------------------------------------------------------------------------------------------------------------"""
    seqlen = len(xios.sequence)
    cut = seqlen // 2
    lines = [[0, seqlen]]
    select = []

    while lines:
        left, right = lines.pop()
        cut = (right - left) // 2
        lines.append([left, cut])
        lines.append([cut + 1, right])

        first = last = None
        for n, s in enumerate(xios.stem_list):
            cross = s.lbegin <= cut and s.rend >= cut
            if cross and first:
                last = n
            elif cross:
                first = n
                last = n

        nstem = last - first + 1
        if nstem <= 25 and nstem > 3:
            # good size selection, this part is done

            n += 1

    return


def overlap(xios):
    """-----------------------------------------------------------------------------------------------------------------
    find the stems that cross each coordinate position

    :param xios: XIOS     RNA structure
    :return:
    -----------------------------------------------------------------------------------------------------------------"""
    pos = [[] for _ in range(len(xios.sequence))]
    stemn = 0
    for stem in xios.stem_list:
        for p in range(stem.lbegin, stem.rend + 1):
            pos[p].append(stemn)
        stemn += 1

    filtered = []
    current = None
    for i, p in enumerate(pos):
        if p != current:
            filtered.append([i, p])
            current = p

    return pos


def stem_overlap(xios):
    """-----------------------------------------------------------------------------------------------------------------
    condensed adjacency matrix. The set of stems that have a non-S relationship to each stem

    :param xios: XIOS     RNA structure
    :return:
    -----------------------------------------------------------------------------------------------------------------"""
    adj = xios.adjacency
    condensed = []
    for i in range(len(adj)):
        condensed.append(set())
        for j in range(len(adj[i])):
            # if adj[i][j] in 's-':
            if adj[i][j] == 's':
                continue
            condensed[-1].add(j)

    return condensed


def common(vselect, vset):
    """
    vset is condensed adjacency matrix - the set of all stems overalapping each stem

    :param vselect:
    :param vset:
    :return:
    """
    common = vset[0]
    for v in list(vselect):
        common = common.intersection(vset[v])
        if not common:
            # no common vertices
            return common

    return common


def cut(piece, branches):
    """

    :param piece:
    :param branches:
    :return:
    """
    parts = []
    incommon = common(piece, branches)
    sbranch = sorted(branches, key=lambda x: len(x), reverse=True)
    used = set()
    for v in sbranch:
        if v.intersection(used):
            # skip if present in used
            continue
        parts.append(v)
        used = used.union(v)
        if used == piece:
            break

    return parts


def extend(v, branches):
    """-----------------------------------------------------------------------------------------------------------------
    extend a vertex set, the extended set is the union of all of its member branches
    :param vset: 
    :param branches: 
    :return: 
    -----------------------------------------------------------------------------------------------------------------"""
    merged = branches[v]
    for thisv in list(branches[v]):
        if thisv > v:
            merged = merged.union(branches[v])
        else:
            merged.add(thisv)

    return merged


def merge(segment, limit):
    """-----------------------------------------------------------------------------------------------------------------
    compare segments and combine when the segments intersect and the union of the two segments does not exceed the
    size limit. Greedy approach: starting from the smallest set, add to all larger sets that have an intersection,
    subject to the size limit.

    :param segment:list of set      vertex subsets from segmentation
    :param limit: int               maximum number of elements in segments sets
    :return: list                   list of segment sets
    -----------------------------------------------------------------------------------------------------------------"""
    final = []
    current = sorted(segment, key=lambda x: len(x), reverse=True)
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


def segment_branch(xios, max):
    """-----------------------------------------------------------------------------------------------------------------

    sort vertices by subtree size
    add all vertices to available
    while available:
        choose largest subtree
        extend
        if less than cutoff
            add to segment
            remove from available

    :param xios: XIOS     RNA structure
    :param max: 
    :return: 
    -----------------------------------------------------------------------------------------------------------------"""
    segment = []
    available = set([i for i in range(len(xios.adjacency))])
    branches = stem_overlap(xios)
    order = [i for i in range(len(branches))]
    sbranch = sorted(branches, key=lambda x: len(x))
    sorder = sorted(order, key=lambda x: len(branches[x]))
    while available:
        try:
            bigv = sorder.pop()
            biggest = branches[bigv]
        except IndexError:
            break
        if len(biggest) > max:
            continue

        if not biggest.intersection(available):
            # there are no unique vertices here
            continue

        extended = extend(bigv, branches)
        if len(extended) > max:
            continue

        segment.append(extended)
        # remove from available
        available.difference_update(extended)
        print(available)

    segment = merge(segment, max)
    return segment


def segment_chunk(xios):
    """-----------------------------------------------------------------------------------------------------------------
    select stems that cross a chunk of the sequence

    :param xios: XIOS     RNA structure
    :return:
    -----------------------------------------------------------------------------------------------------------------"""
    seqlen = len(xios.sequence)
    cut = seqlen // 2
    lines = [[0, seqlen]]
    candidate = []
    condensed = stem_overlap(xios)
    pos = overlap(xios)

    while lines:
        [left, right] = lines.pop()
        half = (right + left) // 2

        first = last = None
        for n, s in enumerate(xios.stem_list):
            l_inchunk = s.lbegin >= left and s.lbegin <= right
            r_inchunk = s.rend >= left and s.rend <= right
            cross = l_inchunk or r_inchunk

            if cross and first:
                last = n
            elif cross:
                first = n
                last = n
        if not first: continue

        nstem = last - first + 1
        if nstem > 25:
            # too big, subdivide
            lines.append([left, half])
            lines.append([half + 1, right])
            continue
        elif nstem <= 3:
            # too small, drop segment
            continue

        candidate.append([first, last])

    # expand these ranges
    adj = xios.adjacency
    for s in candidate:
        first, last = s
        blocklen = last - first + 1
        vlist = set(range(first, last + 1))
        for j in range(0, first):
            aa = adj[j][first:last + 1]
            scount = aa.count('s')
            if scount < blocklen:
                vlist.add(j)
        for j in range(last + 1, len(xios.stem_list)):
            aa = adj[j][first:last + 1]
            scount = aa.count('s')
            if scount < blocklen:
                vlist.add(j)
        print(f'candidate {first},{last} finished \t{len(vlist)}:{vlist}')

    return


# ======================================================================================================================
# Main
# ======================================================================================================================
if __name__ == '__main__':
    sourcedir = sys.argv[1] + '/*.xios'

    for source_xios in glob.glob(sourcedir):
        nfail = 0
        # source_xios = '../data/curated_xios/16S_e.Balamuthia_mandrillaris.xios'
        # print(f'RNA:{source_xios}')
        rna = Topology(source_xios, xios=source_xios)
        segment = segment_branch(rna, 25)
        for s in segment:
            if len(s) > 25:
                print(f'{source_xios} failed')
                nfail = 0

    print(f'number failed:{nfail}')

    exit(0)
