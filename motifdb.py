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

    def addstemall(self):
        """-----------------------------------------------------------------------------------------
        return the set of structures with one additional stem explicitly added at all possible
        positions. This is just to check the other enumeration approaches.

        :return: list of SerialRNA
        -----------------------------------------------------------------------------------------"""
        # make a new structure with the stem number incremented by 1
        newlen = len(self) + 2
        children = []

        for begin in range(0, newlen - 1):
            for end in range(begin + 1, newlen):
                extended_rna = [None for _ in range(newlen)]
                extended_rna[begin] = 0
                extended_rna[end] = 0
                newpos = 0
                for pos in self:
                    while extended_rna[newpos] is not None:
                        newpos += 1
                    extended_rna[newpos] = pos + 1

                children.append(SerialRNA(extended_rna))

        return children

    def addstemleft(self):
        """-----------------------------------------------------------------------------------------
        return the set of structures with one additional stem added at all possible position

        :return: list of SerialRNA
        -----------------------------------------------------------------------------------------"""
        # make a new structure with the stem number incremented by 1
        base = []
        newlen = len(self) + 2
        half = newlen // 2

        children = []
        for pos in self:
            base.append(pos + 1)

        for begin in range(0, half + 1):
            for end in range(begin + 1, newlen):
                extended_rna = [None for _ in range(len(self) + 2)]
                extended_rna[begin] = 0
                extended_rna[end] = 0
                newpos = 0
                for pos in base:
                    while extended_rna[newpos] is not None:
                        newpos += 1
                    extended_rna[newpos] = pos

                children.append(SerialRNA(extended_rna))

        return children

    def addstemzero(self):
        """-----------------------------------------------------------------------------------------
        return the set of structures with one additional stem such that the beginning of the stem
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

        :return: True
        -----------------------------------------------------------------------------------------"""
        stem = 0
        map = {}
        for pos in self:
            if pos not in map:
                map[pos] = stem
                stem += 1

        for pos in range(len(self)):
            self[pos] = map[self[pos]]

        return True

    def reverse(self):
        """-----------------------------------------------------------------------------------------
        reverse the rna structure, numbering is not converted to canonical
        :return: True
        -----------------------------------------------------------------------------------------"""
        begin = 0
        end = len(self) - 1
        while begin < end:
            self[begin], self[end] = self[end], self[begin]
            begin += 1
            end -= 1

        return True

    def canonical_fbstr(self):
        """-----------------------------------------------------------------------------------------
        return the graph and its reverse as strings in canonical form

        :return: str, str, forward and backward canonical strings
        -----------------------------------------------------------------------------------------"""
        self.canonical()
        forward = self.tostring()

        self.reverse()

        self.canonical()
        backward = self.tostring()

        return forward, backward

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
            #typically the reversed structure is absent
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
    print('elapsed time {:.4f}'.format(stop-start))

    exit(0)
