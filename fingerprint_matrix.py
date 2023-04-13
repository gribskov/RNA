"""=================================================================================================
Convert a set of fingerprint files into a binary form for rapid reading and manipulation
also does a kmeans clustering based on the motifs present in each fingerprint and identifies the
most correlateed motifs

usage:
    python fingerprint_matrix.py <input fileglob>
================================================================================================="""
import sys
from math import floor, ceil

from fingerprint import FingerprintMatrix
from kmeans import Kmeans

####################################################################################################
# main program
####################################################################################################
if __name__ == '__main__':
    n_repeat = 150
    k = Kmeans(30)
    neighborhood_cutoff = (0.4, 0.6)
    show_cycle = False
    show_init = False

    selection = sys.argv[1]
    print( f'selected files: {selection}')
    fmat = FingerprintMatrix()
    fmat.read_files(selection)

    fmatfile = 'fmatrix.tsv')
    fmat.write('fmatfile'
    print(f'fingerprint matrix written to {fmatfile}')

    fmatpkl = 'fmatrix.pl'
    fmat.pickle(fmatpkl)
    print(f'fingerprint matrix written to {fmatpkl}')

    # fmat.select_min_max(10, 45, False, recalculate=True)

    together = [[0 for _ in range(len(fmat.fpt))] for _ in range(len(fmat.fpt))]

    for repeat in range(n_repeat):
        # cluster n_repeat times
        # TODO save the lowest error clustering
        ndata = k.assign_data_random(fmat.fpt)
        if show_init:
            print(f'{ndata} points from {selection} for k={k.k}')
            print(f'initial groups')
            for g in range(len(k.group)):
                print(f'group {g}: ', end='')
                for i in k.group[g]:
                    print(f'\t{i}', end='')
                print()

        print(f'repeat {repeat}')
        k.cluster()

        g = 0
        fptname = list(fmat.fpt_id.keys())
        for group in k.group:
            g += 1
            if show_cycle:
                # print the clusters at each cycle
                print(f'group {g}')
                for f in group:
                    print(f'\t{fptname[f]}')

            pairs = ((a, b) for idx, a in enumerate(group) for b in group[idx + 1:])
            for a, b in pairs:
                together[a][b] += 1
                together[b][a] += 1

        # final clusters
        # g=0
        # for group in k.group:
        #     g += 1
        #     print(f'group {g}')
        #     for f in group:
        #         print(f'\t{fptname[f]}')

    # neighborhood clustering: group the points that co-occur the most often using single linkage
    # clustering
    # for minval in range(floor(n_repeat * neighborhood_cutoff), n_repeat):
    for minval in range(floor(n_repeat * neighborhood_cutoff[0]),
                        ceil(n_repeat * neighborhood_cutoff[1])):
        print(f'\n{"#" * 80}\n minval={minval}\n{"#" * 80}')
        cluster = []
        n = 0
        avail = [i for i in range(len(fmat.fpt_id))]
        while avail:
            c = avail.pop()
            cluster.append([c])
            todo = [c]
            while todo:
                p0 = todo.pop()
                for p1 in avail:
                    if together[p0][p1] > minval:
                        cluster[n].append(p1)
                        todo.append(p1)
                        avail.remove(p1)

            n += 1

        c = 0
        for neighborhood in cluster:
            print(f'cluster {c}')
            c += 1
            for f in neighborhood:
                print(f'\t{fptname[f]}')

exit(0)
