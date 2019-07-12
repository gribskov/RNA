"""=================================================================================================
sample a random fingerprint from a structure in XIOS XML format

10 July 2019     Michael Gribskov
================================================================================================="""
import sys
from topology import Topology
from xios import Xios, Gspan, MotifDB
from fingerprint import Fingerprint

# read in the database
motif_file = 'data/5stem.list.txt'
motifdb = MotifDB(json=motif_file)

# read in the structure as XIOS XML
structurefile = 'data/rnasep_a1.Buchnera_APS.xios'
rna = Topology(xml=structurefile)

fingerprint = Fingerprint()
subgraph = 5

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

print('Simple fingerprint: {}\t{}\t{}'.format(fingerprint.total, fingerprint.nmotif, fingerprint.mincount()))
fingerprint.add_parents(motifdb)
print('Extended fingerprint: {}\t{}\t{}'.format(fingerprint.total, fingerprint.nmotif, fingerprint.mincount()))

print( fingerprint.toYAML())
fingerprint.writeYAML('data/test.fpt')

exit(0)
