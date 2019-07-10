"""=================================================================================================
sample a random fingerprint from a structure in XIOS XML format

10 July 2019     Michael Gribskov
================================================================================================="""
from topology import Topology
from xios import Xios, Gspan

rna = Topology()
rna.XIOSread('data/rnasep_a1.Buchnera_APS.xios')

fingerprint = {}
nrandom = 0
nmotif = 0
block = 1000
# for i in range(50000):
mincount = 0
while mincount < 3:
    for i in range(block):
        nrandom += 1
        xios = rna.sample_xios(5)
        gspan = Gspan(graph=xios)
        dfs = gspan.minDFS().ascii_encode()

        if dfs in fingerprint:
            fingerprint[dfs] += 1
        else:
            fingerprint[dfs] = 1
            nmotif += 1

    mincount = min(fingerprint.values())
    print('{}\t{}\t{}'.format(nrandom, nmotif, mincount))

n = 0
for dfs in sorted(fingerprint, key=lambda dfs: fingerprint[dfs], reverse=True):
    n += 1
    print('{:5}\t{:6}\t{:7.5f}\t{}'.format(n, fingerprint[dfs], fingerprint[dfs]/nrandom, dfs))

exit(0)