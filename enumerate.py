"""=================================================================================================
Enumerate all possible XIOS motifs

biological structures are enumerated by adding new edges to previousl enumerated structures at all
possible positions.  Many of the structures generated in this way are duplicates or inversions of
previously processed structures (inverted structures should have the same minimum DFS code).

Points to consider:
Unique minimum DFS identifies the unique motifs, but is expensive to calculate using gspan

The string form of the the SerialRNA form can be used to detect many duplicates and inversions, but
some isomorphs can only be identified with gspan.

Unconnected structures can be joined together by wrapping them in a single new stem.  This means
that all structures, connected and unconnected, with n stems must be saved for generating
structures withn+1 stems.

work with three lists:
current: structures that can be extended to generate larger strucutres (includes unconnected)
allstr: string format of canonically numbered structures that have been processed.  Used to identify
        structures that have alreayd been processed.  Both forward and backward orientations are
        included.  This is implemented as a hash to keep only unique strings
motif: structures with unique minimum DFS codes

push starting graph on current

while current
    pop rna from beginning of current
    generate extensions of rna -> candidate list
    remove rna from current

    for each rna in candidate
        check forward and backward canonical SerialRNA strings against allstr list
        if  fstr or bstr present, break

        add fstr and bstr to allstr
        add rna to end of current
        if rna is connected
            find minimum DFS code with gspan
            if unique add to motifs

================================================================================================="""
import time
from motifdb import SerialRNA
from xios import Xios
from gspan import Gspan

maxgraph = 12
current = [SerialRNA([0, 0])]
candidate = []
motif = {}
allstr = {}

start = time.time()
while True:
    # pop a structure from the list of RNAs to be extended
    thisrna = current[0]
    # print('rna {}     candidate {}'.format(thisrna, candidate))

    if len(thisrna) >= maxgraph:
        # print(current)
        break

    # generate extended structures, and remove current structure from list of RNAs to be extended
    candidate = thisrna.addstemleft()
    current.remove(thisrna)

    for rna in candidate:

        fstr, bstr = rna.canonical_fbstr()
        if fstr in allstr or bstr in allstr:
            # skip this structure if it has been seen before
            continue

        # update list of previously processed structures
        allstr[fstr] = 1
        allstr[bstr] = 1

        # save both unconnected and unique connected for the next round of extension
        current.append(rna)

        if len(rna.connected()) == 1:
            xios = Xios(serial=rna)
            gspan = Gspan(xios)
            dfs = gspan.minDFS()
            dfshex = dfs.ascii_encode()
            if dfshex not in motif:
                # save unique minimum DFS codes
                motif[dfshex] = {'str': fstr, 'min': dfs}
                # print('\tmotif {}\t{}'.format(fstr, motif[dfshex]))

# o = open('data/12stem.list.txt', 'w')
stop = time.time()
print('motifs {}'.format(len(motif)))
for m in motif:
    # if len(motif[m]['str'])>12:
    #     o.write('{}\t{}\t{}\n'.format(m, motif[m]['str'], motif[m]['min']))
    print('{}\t{}\t{}'.format(m, motif[m]['str'], motif[m]['min']))
print('motifs {}'.format(len(motif)))
print('elapsed time {:.4f}'.format(stop - start))

exit(0)
