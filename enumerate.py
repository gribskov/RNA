"""=================================================================================================
Enumerate all possible XIOS motifs

================================================================================================="""
import time
from xios import Xios
from graph import RNAGraph
from gspan import Gspan

rna = SerialRNA([0, 0])
graph = RNAGraph(rna)
print(graph.toList())
gg = graph.toList()
xios = Xios()
xios.from_graph(graph.pairs)
gspan = Gspan(xios)

# unique = {'string': rna.tostring(), 'structure': rna, 'mindfs': gspan.minDFS()}
current = [rna]
candidate = []
motif = {}
allstr = {}
maxgraph = 14

start = time.time()
while True:
    thisrna = current[0]
    if len(thisrna) >= maxgraph:
        # print(current)
        break

    candidate = thisrna.addstemall()
    try:
        current.remove(thisrna)
    except ValueError:
        # typically the reversed structure is absent
        pass

    # print('rna {}     candidate {}'.format(thisrna, candidate))
    for rna in candidate:

        # print('{}'.format(thisrna))
        # rna.canonical()
        # graphstr = rna.tostring()
        # if graphstr in allstr:
        #     # print('skipping {}'.format(graphstr))
        #     continue
        #
        # allstr.append(graphstr)
        fstr, bstr = rna.canonical_fbstr()
        if fstr in allstr or bstr in allstr:
            continue

        allstr[fstr] = 1
        allstr[bstr] = 1

        graphstr = fstr
        if len(rna.connected()) == 1:
            graph = RNAGraph(rna)
            xios = Xios()
            xios.from_graph(graph.pairs)
            gspan = Gspan(xios)
            dfs = gspan.minDFS()
            dfsxios = Xios()
            dfsxios.from_list(dfs)
            dfshex = dfsxios.ascii_encode()
            if dfshex not in motif:
                motif[dfshex] = {'str': graphstr, 'min': dfs}
                # print('\tmotif {}\t{}'.format(graphstr, motif[dfshex]))

            else:
                # this is a duplicate, remove from candidate, do not save on current
                # print('\tduplicate {}'.format(rna))
                # candidate.remove(rna)
                continue

        else:
            # print('\tnot connected')
            pass

        # save unconnected and unique connected
        current.append(rna)

stop = time.time()
print('motifs {}'.format(len(motif)))
for m in motif:
    print('{}\t{}\t{}'.format(m, motif[m]['str'], motif[m]['min']))
print('motifs {}'.format(len(motif)))
print('elapsed time {:.4f}'.format(stop - start))

exit(0)
