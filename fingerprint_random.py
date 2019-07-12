"""=================================================================================================
sample a random fingerprint from a structure in XIOS XML format

10 July 2019     Michael Gribskov
================================================================================================="""

# other sampling strategies, tried but didn't work better
# limit = 1
# last = 0
# scale = 1.5
# gap = 1
# while nrandom < limit:
#     nrandom += 1
#     xios = rna.sample_xios(5)
#     gspan = Gspan(graph=xios)
#     dfs = gspan.minDFS().ascii_encode()
#
#     if dfs in fingerprint:
#         fingerprint[dfs] += 1
#     else:
#         fingerprint[dfs] = 1
#         nmotif += 1
#         gap = max( gap, nrandom-last)
#         limit = nrandom + gap * scale
#
#     mincount = min(fingerprint.values())
#     if mincount > 3:
#         break
#
#     print('{}\t{}\t{}\t{}\t{}'.format(nrandom, nmotif, mincount, gap, limit))

# set all waits to intial weight, p (a +1 prior)
# w = []
# recip = []
# p = 1
# for stem in range(len(rna.stem_list)):
#     w.append(p)
#     recip.append(0)
#
# mincount = 0
# while mincount < 2:
#     for i in range(block):
#         # reciprocal weights
#         # for i in range(len(w)):
#             # recip[i] = 1.0 / w[i]
#         # order =  sorted([i for i in range(len(w))], key=lambda x:w[x], reverse=True)
#         r = 1
#         for idx in sorted([i for i in range(len(w))], key=lambda x:w[x], reverse=True):
#             recip[idx] = r
#             r += .25
#         nrandom += 1
#         # print(recip)
#         # print(w, '\n')
#         xios, vlist = rna.sample_xios_weighted(5, recip)
#         for v in vlist:
#             w[v] += 1
#         gspan = Gspan(graph=xios)
#         dfs = gspan.minDFS().ascii_encode()
#
#         if dfs in fingerprint:
#             fingerprint[dfs] += 1
#         else:
#             fingerprint[dfs] = 1
#             nmotif += 1
#
#     mincount = min(fingerprint.values())
#     print('{}\t{}\t{}'.format(nrandom, nmotif, mincount))

# # weight by number of neighbors
# adj = rna.adjacency
# ncount = []
# for row in range(len(adj)):
#     ncount.append(0)
#     for col in range(len(adj)):
#         if adj[row][col] in ('i', 'j', 'o' ):
#             ncount[row] += 1
#
#
#
# mincount = 0
# exp = 1.0
# w = [ 0 for _ in range(len(ncount))]
# while mincount < 2:
#     print(ncount)
#     for i in range(len(ncount)):
#         w[i] = ncount[i] ** exp
#     print(w)
#     exp -= 0.1
#     for i in range(block):
#
#         nrandom += 1
#         # print(recip)
#         # print(w, '\n')
#         xios, vlist = rna.sample_xios_weighted(6, w)
#         gspan = Gspan(graph=xios)
#         dfs = gspan.minDFS().ascii_encode()
#
#         if dfs in fingerprint:
#             fingerprint[dfs] += 1
#         else:
#             fingerprint[dfs] = 1
#             nmotif += 1

import sys
from topology import Topology
from xios import Xios, Gspan
from fingerprint import Fingerprint

# read in the database
motiffile = 'data/5stem.list.txt'
try:
    motif = open(motiffile, 'r')
except OSError:
    sys.stderr.write('fingerprint_random - unable to open motif dictionary ({}'.format(motiffile))

# read in the structure as XIOS XML
structurefile = 'data/rnasep_a1.Buchnera_APS.xios'
rna = Topology(xml=structurefile)

fingerprint = Fingerprint()
subgraph = 5
flist = []
nrandom = 0
block = 1000

# unweighted seems to be the best
max = []
cycles = []
count_threshold = 2

xios = rna.sample_xios(subgraph)
gspan = Gspan(graph=xios)
dfs = gspan.minDFS().human_encode()
fingerprint.add(dfs)
minmotif = fingerprint.minkey()

while True:
    # for i in range(block):
    xios = rna.sample_xios(subgraph)
    gspan = Gspan(graph=xios)
    dfs = gspan.minDFS().human_encode()
    fingerprint.add(dfs)

    if dfs == minmotif:
        minmotif = fingerprint.minkey()
        mincount = fingerprint.mincount()
        # print('{}\t{}\t{}'.format(fingerprint.total, fingerprint.nmotif, fingerprint.mincount()))
        if mincount > count_threshold:
            break

print('{}\t{}\t{}'.format(fingerprint.total, fingerprint.nmotif, fingerprint.mincount()))
n = 0
# for dfs in sorted(fingerprint, key=lambda dfs: fingerprint[dfs], reverse=True):
#     n += 1
#     print('{:5}\t{:6}\t{:9.5f}     {}'.format(n, fingerprint[dfs], fingerprint[dfs]/nrandom, dfs))

exit(0)
