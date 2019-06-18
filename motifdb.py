"""=================================================================================================
MotifDB class is for creating and using XIOS graph dictionaries

011220

0,5  1,2  3,4

00 -> 11 *
    0011 -> 1122
        010122
        011022
        011202 *
        011220 *
    0101 -> 1212  *
        010212 *
        012012 *
        012102 *
        012120 *
    0110 -> 1221  *
        010221 *
        012021 *
        012201 *
        012210 *

01
    01 23

    02 13
    03 12


Michael Gribskov     15 June 2019
================================================================================================="""


class MotifDB(list):
    """=============================================================================================

    ============================================================================================="""

    # __init__ inherited from list
    def dummy(self):
        pass
        return


class SerialRNA(list):
    """=============================================================================================
    for working with RNAS encoded for example as 001212, meaning ( ) ( [ ) ]
    ============================================================================================="""

    def connected(self):
        """-----------------------------------------------------------------------------------------
        return the connected graph(s) in the current structure. The entire graph is connected if one
        structure is returned

        :return: list of SerialRNA
        -----------------------------------------------------------------------------------------"""
        component = []
        open = []
        begin = 0
        for pos in range(len(self)):
            if self[pos] in open:
                open.remove(self[pos])
                if len(open) == 0:
                    component.append(SerialRNA(self[begin:pos + 1]))
                    begin = pos + 1
            else:
                open.append(self[pos])

        if len(component) == 1:
            return [self]
        else:
            return component

    def addstemleft(self):
        """-----------------------------------------------------------------------------------------
        return the set of structures with one additional stem added at all possible position

        :return: list of SerialRNA
        -----------------------------------------------------------------------------------------"""
        # make a new structure with the stem number incremented by 1
        base = []
        children = []
        for pos in self:
            base.append(pos + 1)

        for end in range(0, len(base)):
            extended_rna = [0]
            for pos in range(end):
                extended_rna.append(base[pos])
            extended_rna.append(0)
            for pos in range(end, len(base)):
                extended_rna.append(base[pos])
            children.append(SerialRNA(extended_rna))

        return children

    def addstemzero(self):
        """-----------------------------------------------------------------------------------------
        return the set of structures with one additional stem duvh that the beginning of the stem
        is at position zero

        :return: list of SerialRNA
        -----------------------------------------------------------------------------------------"""
        # make a new structure with the stem number incremented by 1
        base = []
        children = []
        newlen = len(self) + 2
        for pos in self:
            base.append(pos + 1)

        for end in range(1, newlen):
            extended_rna = [None for _ in range(newlen)]
            extended_rna[0] = 0
            extended_rna[end] = 0
            basepos = 0
            for pos in range(1, newlen):
                if extended_rna[pos] is None:
                    extended_rna[pos] = base[basepos]
                    basepos += 1

            children.append(SerialRNA(extended_rna))

        return children

    def canonical(self):
        """-----------------------------------------------------------------------------------------
        convert graph to canonical form.  In canonical form the stems occur in increasing numerical
        order beginning at zero

        :return: True if graph is changed, otherwise False
        -----------------------------------------------------------------------------------------"""
        stem = 0
        changed = False
        map = {}
        for pos in self:
            if pos not in map:
                map[pos] = stem
                if pos != stem:
                    changed = True
                stem += 1
        if changed:
            for pos in range(len(self)):
                self[pos] = map[self[pos]]

        return changed

    def tostring(self):
        """-----------------------------------------------------------------------------------------
        return a string representing the stucture.  the string is the concatenation of the digits.

        :return: string
        -----------------------------------------------------------------------------------------"""
        return ''.join(str(x) for x in self)


# --------------------------------------------------------------------------------------------------
# testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    rnas = [[0, 0, 1, 1, 2, 2],
            [0, 1, 0, 1, 2, 2],
            [0, 1, 1, 2, 2, 0],
            [0, 1, 2, 1, 2, 0],
            [0, 0], []
            ]

    print('canonical form')
    noncanonical = [[3, 3, 0, 0, 1, 1], [1, 1, 2, 2, 3, 3], [1, 1, 2, 2, 4, 4], [3, 2, 0, 2, 0, 3]]
    for testcase in noncanonical:
        rna = SerialRNA(testcase)
        print('RNA {}'.format(rna))
        rna.canonical()
        print('\tcanonical {}'.format(rna))

    print('\nConnected components')
    for testcase in rnas:
        rna = SerialRNA(testcase)
        print('RNA {}'.format(rna))
        connected = rna.connected()
        if len(connected) > 1:
            for i in range(len(connected)):
                print('\tcomponent {}: {}'.format(i, connected[i]))

    print('\nExtension')
    for testcase in rnas:
        rna = SerialRNA(testcase)
        print('RNA {}'.format(rna))
        for new in rna.addstemzero():
            print('\t{} {}'.format(new, len(new.connected())))

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
    maxgraph = 14

    while True:
        thisrna = current[0]
        if len(thisrna) >= maxgraph:
            break

        candidate = thisrna.addstemzero()
        current.remove(thisrna)

        # print('\nrna {}     candidate {}'.format(thisrna, candidate))

        for thisrna in candidate:
            print('{}'.format(thisrna))
            graphstr = thisrna.tostring()
            if len(thisrna.connected()) == 1:
                # print('\tconnected')
                graph = RNAGraph(thisrna)
                xios = Xios()
                xios.from_graph(graph.pairs)
                gspan = Gspan(xios)
                dfs = gspan.minDFS()
                dfsxios = Xios()
                dfsxios.from_list(dfs)
                dfshex = dfsxios.hex2_encode()
                if dfshex not in motif:
                    motif[dfshex] = {'str': graphstr, 'min': dfs}
                    # print('\tmotif {}\t{}'.format(graphstr, motif[dfshex]))

                else:
                    # this is a duplicate, remove from candidate, do not save on current
                    # print('\tduplicate')
                    candidate.remove(thisrna)
                    continue
            else:
                # print('\tnot connected')
                pass

            # save unconnected and unique connected
            current.append(thisrna)

    print( 'motifs {}'.format(len(motif)))
    for m in motif:
        print('{}\t{}\t{}'.format(m, motif[m]['str'], motif[m]['min']))
    print('motifs {}'.format(len(motif)))

    exit(0)
