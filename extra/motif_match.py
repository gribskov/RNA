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
        for v in range(len(adj)):
            e = adj[v][u]
            if e not in 'ijo':
                continue

            idx[e].append([u,v])

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

    map = [[[None for _ in range(vmax+1)],0]]
    while map:
        # depth first search of d to graph mappings
        mapping, edge_n = map.pop()
        if edge_n == len(motif):
            print(mapping)
            continue
        d0, d1, edge = motif[edge_n]
        for g0, g1 in idx[edge]:
            # all possible matching edges to the DFS row

            # check for match already defined vertices
            if mapping[d0] is not None:
                # vertex is known
                if g0 != mapping[d0]:
                    # discordant with current mapping
                    continue

            try:
                if mapping[d1] is not None:
                    # vertex is known
                    if g1 != mapping[d1]:
                        # discordant with current mapping
                        continue
            except IndexError:
                print('oops')

            # try:
            #     g0idx = mapping.index(g0)
            # except ValueError:
            #     #g0 not found in map
            #     g0idx = None
            # if g0idx != d0 and g0idx is not None:
            #     continue
            #
            # try:
            #     g1idx = mapping.index(g1)
            # except ValueError:
            #     #g0 not found in map
            #     g1idx = None
            # if g1idx != d0 and g1idx is not None:
            #     continue



            # if g0 in mapping or g1 in mapping:
            #     # g0 and g1 are already in current mapping
            #     continue
            #
            # if mapping[d0] != None and mapping [d0] != g0:
            #     # d0 and d1 are already defined
            #     continue
            # if mapping[d1] != None and mapping [d1] != g1:
            #     # d0 and d1 are already defined
            #     continue

            new = [mapping[:],edge_n+1]
            new[0][d0] = g0
            new[0][d1] = g1
            map.append(new)

    exit(0)
