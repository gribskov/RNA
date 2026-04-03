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
            segment.append([begin,end])
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
        position.append([s.lbegin,s.rend, n, 1])
        position.append([s.rend, s.lbegin, n, -1])
        n += 1

    height = 0
    n = 0
    for pos in sorted(position, key=lambda x:x[0]):
        height += pos[3]
        print(f'{n}\t{height}\t{pos}')
        n+= 1

    return

# ======================================================================================================================
# Main
# ======================================================================================================================
if __name__ == '__main__':
    sourcedir = sys.argv[1] + '/*.xios'

    for source_xios in glob.glob(sourcedir):
        print(f'{source_xios}')
        rna = Topology(source_xios,xios=source_xios)
        segment = segment_mountain(rna)
        print(f'\t{source_xios}\t{len(rna.stem_list)}\t{len(segment)}\t{segment}')

    exit(0)
