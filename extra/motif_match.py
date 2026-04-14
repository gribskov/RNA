"""=====================================================================================================================
motif_match.py

find matches of a motif in a xios structure
motif example: 0i1.1i2.2j0.1o3.3j0.0i4.0i5.0i6.

Michael Gribskov 4/14/2026
====================================================================================================================="""
import sys
from xios import Xios
from topology import Topology

def edge_parse(motifstr):
    """---------------------------------------------------------------------------------------------
    convert motif string to a list of edges. each edge is [v0, v1, edge]

    :param motif:
    :return:
    ---------------------------------------------------------------------------------------------"""
    motif = motifstr.split('.')
    motif.remove('')
    parsed = []
    big_v = 0
    for m in motif:
        # find max vertex
        if 'i' in m:
            e = 'i'
        elif 'j' in m:
            e = 'j'
        elif 'o' in m:
            e = 'o'
        else:
            sys.stderr.write(f'unknown edge in {m}')

        # print(m)
        v0, v1 = m.split(e)
        v0 = int(v0)
        v1 = int(v1)
        big_v = max(big_v, v0, v1)
        parsed.append([v0, v1, e])

    return big_v, parsed

def edge_index(adj):
    """---------------------------------------------------------------------------------------------
    make index of i, j, o edges values are vertex numbers

    :param adj:
    :return:
    ---------------------------------------------------------------------------------------------"""
    idx = {'i':[], 'j':[], 'o':[]}

    for u in range(len(adj)):
        for x in 'ijo':
            idx[x].append([])
        for v in range(len(adj)):
            e = adj[u][v]
            if e not in 'ijo':
                continue

            idx[e][u].append(v)

    return idx


# ======================================================================================================================
# main
# ======================================================================================================================
if __name__ == '__main__':
    xios = Topology(xml=sys.argv[1])
    motifstr = sys.argv[2]

    vmax, motif = edge_parse(motifstr)
    d2g = [[] for _ in range(vmax+1)]

    idx = edge_index(xios.adjacency)

    # push initial mappings based on first DFS row onto map stack
    d0, d1, edge = motif[0]
    map = []
    for u,vlist in enumerate(idx[edge]):
        for v in vlist:
            mapping = [None for _ in range(vmax+1)]
            mapping[d0] = u
            mapping[d1] = v
            map.append([mapping,1])

    while map:
        # depth first search of d to graph mappings
        mapping, edge_n = map.pop()
        if edge_n == len(motif):
            print(mapping)
            continue
        d0, d1, edge = motif[edge_n]
        g0 = mapping[d0]
        if mapping[d1] is not None:
            # both d0 and d1 are defined, check if the required edge exists
            # this is a backward edge
            if mapping[d1] in idx[edge][g0]:
                map.append([mapping,edge_n + 1])
        else:
            # d1 is not yet define, push mappings for all matches
            # forward edges
            # print(idx[edge][g0])
            for v in idx[edge][g0]:
                if v not in mapping:
                    new = mapping[:]
                    new[d1] = v
                    map.append([new,edge_n + 1])


    exit(0)
