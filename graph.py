class RNAGraph:
    """=============================================================================================
    RNA graph class

    There are several ways to represent a folded RNA structure as a linear string at an abstract
    level where we consider simply paired regions rather than individual base pairs.  These include
    the Giegerich abstract shapes approach, where the runs of parentheses and dots in a Vienna
    formatted structure are condensed to single brackets (square brackets in their formulation).
    This package handles the manipulation of these kinds of linear representations (including some
    extensions for pseudoknots). A second representation numbers the stem regions based on the
    pairing.  Again, this can be condensed so each stem is represented by a single digit
    For instance, for a tRNA

    >S.cerevisiae_tRNA-PHE M10740/1-73
    GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA
    (((((((..((((........)))).((((.........)))).....(((((.......)))))))))))).   Vienna format
    1111111  2222        2222 3333         3333     44444       444441111111    stem number
    1111111  2222        3333 4444         5555      6666       777788888888    position

    abstract shape (level 3) [ [ ] [ ] [ ] ]    (spaces added for readability)
    list format   1 2 2 3 3 4 4 1               index indicates position, value indicates stem
    pair format   (1,8) (2,3) (4,5) (6,7)       each item is a stem, values are the positions
                  1 8 2 3 4 5 6 7               (serialized pair format)

    example with a pseudoknot added
    GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA
    (((((((..((((........)))).((((.........)))).[[[.(((((..]]]..)))))))))))).   Vienna format
    1111111  2222        2222 3333         3333 444 55555  444  555551111111    stem number
    1111111  2222        3333 4444         5555 666 77777  888  999991111111    position
                                                                     0000000
    ( ( ) ( ) [ ( ] ) )
    1 2 2 3 3 4 5 4 5 1
    (1, 10) (2, 3) (4, 5) (6, 8) (7, 9)
    1 10 2 3 4 5 6 8 7 9

    Synopsis
        from graph import RNAGraph

        graph = RNAGraph()
        graph.fromList([0,1,1,0])
            or
        graph = RNAGraph(l=[0,1,1,0])

        graph.reverse()     # reverse the order of stems left-right
    ============================================================================================="""

    def __init__(self, inlist=[]):
        """-----------------------------------------------------------------------------------------
        Graph constructor
        Internal data structure is a pair graph: a list in which each element is a stem.  The stems
        give the coordinates of the left and right half stems as integers along a line.
        -----------------------------------------------------------------------------------------"""
        self.pairs = []
        self.nstem = 0

        if inlist:
            self.fromList(inlist)

    def __str__(self, sep='_'):
        """-----------------------------------------------------------------------------------------
        Return a serialized version of the pair structure
        :return: string
        -----------------------------------------------------------------------------------------"""

        s = ''
        for p in self.pairs:
            s += '{}_{}_'.format(p[0], p[1])

        return s.rstrip('_')

    def __len__(self):
        """-----------------------------------------------------------------------------------------
        length is the number of stems
        -----------------------------------------------------------------------------------------"""
        return len(self.pairs)

    def fromList(self, g):
        """"----------------------------------------------------------------------------------------
        read a graph in list format as a list.  returns a list of lists with the begin/end position
        of each stem (pair format)
        :param g: list, structure in list format
        :return: integer, number of stems
        -----------------------------------------------------------------------------------------"""
        self.nstem = int(len(g) / 2)
        self.pairs = [[] for _ in range(self.nstem)]

        # the values indicate the stem number, the index indicates the position
        for i in range(len(g)):
            self.pairs[g[i]].append(i)

        return self.nstem

    def fromListAsString(self, g, sep=' '):
        """-----------------------------------------------------------------------------------------
        the input list is a string separated by sep
        :param g: string, input graph as a string
        :param sep: string, separation character in input string
        :return: integer, number of stems
        -----------------------------------------------------------------------------------------"""
        listval = g.split(sep)
        nstem = self.fromList(listval)

        return nstem

    def toList(self):
        """-----------------------------------------------------------------------------------------
        return the structure in list format (s a list)
        :return g: list, structure in list format
        -----------------------------------------------------------------------------------------"""
        g = [0 for _ in range(self.nstem * 2)]

        stem = 0
        for pair in self.pairs:
            g[pair[0]] = stem
            g[pair[1]] = stem
            stem += 1

        return g

    def fromVienna(self):
        """-----------------------------------------------------------------------------------------
        Read structure in Vienna format
        :return: integer, number of stems
        -----------------------------------------------------------------------------------------"""
        pass

    def toVienna(self, pad=' '):
        """-----------------------------------------------------------------------------------------
        return a string with the structure in vienna format.  This is basically the abstract shapes
        format with support for pseudoknots.
        :return: string
        -----------------------------------------------------------------------------------------"""
        bracket = [['(', ')'], ['[', ']'], ['{', '}'], ['<', '>'], [':', ':']]
        vienna = ['.' for _ in range(self.nstem * 2)]
        list_format = self.toList()

        level = []
        knot = True
        for stemi in self.pairs:
            avail = 0
            for l in range(len((level))):
                knot = True
                for stem_num in range(len(level[l])):
                    stemj = self.pairs[stem_num]
                    if stemi[0] > stemj[1] or (stemi[0] > stemj[0] and stemi[1] < stemj[1]):
                        # stemi is nested or serial with the stem using this level (stemj)
                        knot = False
                    else:
                        knot = True
                        break

                if not knot:
                    level[l].append(stemi)
                    break
            if knot:
                level.append([])
                level[-1].append(stemi)

        for l in range(len(level)):
            for stem in level[l]:
                vienna[stem[0]] = bracket[l][0]
                vienna[stem[1]] = bracket[l][1]

        return pad.join(vienna)

    def reverse(self):
        """-----------------------------------------------------------------------------------------
        reverse the positions of the stems by converting to maxpos-pos and resorting in order
        of first coordinate
        :return: None
        -----------------------------------------------------------------------------------------"""
        m = self.nstem * 2 - 1
        for i in range(self.nstem):
            self.pairs[i][0], self.pairs[i][1] = m - self.pairs[i][1], m - self.pairs[i][0]

        self.pairs.sort(key=lambda k: k[0])

        return self.pairs

    def connected(self):
        """-----------------------------------------------------------------------------------------
        Returns True if graph is i-o connected
        :return: True / False
        -----------------------------------------------------------------------------------------"""
        pos = 1
        for d in self.depth():
            if pos == 0:
                continue
            pos += 1
            if d == 0:
                break

        if pos < len(self) * 2:
            return False

        return True

    def depth(self):
        """-----------------------------------------------------------------------------------------
        Return  list with the nesting depth at each position of the graph in list format.  This is
        useful for mountain plots and determining connectivity.  If the level reaches zero before
        the last position, the graph is disconnected.
        :return: list
        -----------------------------------------------------------------------------------------"""
        depth = []
        d = 0
        stem = []

        for s in self.toList():
            if s in stem:
                # stem seen before
                d -= 1
            else:
                # new stem
                stem.append(s)
                d += 1

            depth.append(d)

        return depth


# --------------------------------------------------------------------------------------------------
# Non-object functions
# --------------------------------------------------------------------------------------------------
def possible(s):
    """---------------------------------------------------------------------------------------------
    return a list of the possible next stems.  Two kinds of stems can be added
    1) stems with use zero can be added in numerical order
    2) stems that have use == 1 can be added, but only if the the count of open stems would be > 0

    :param s: stem usage list
    :return: list of possible elements to add
    ---------------------------------------------------------------------------------------------"""
    possible = []
    candidate = []
    first = True
    open = 0
    all = 0
    for i in range(len(s)):
        all += 2 - s[i]
        if s[i] == 0 and first:
            first = False
            possible.append(i)
        elif s[i] == 1:
            candidate.append(i)
            open += 1

    if open > 1 or all == 1:
        possible += candidate

    return possible


def enumerateRNATopology(n):
    """---------------------------------------------------------------------------------------------
    enumerate all graphs with n stems
    :param n: number of stems
    :return: array of graph lists
    ---------------------------------------------------------------------------------------------"""
    graphs = []
    s = [0 for _ in range(n)]
    stack = []
    struc = []
    stack.append((0, s[:], struc))

    while stack:
        n, s, struc = stack.pop()
        struc.append(n)
        s[n] += 1

        p = possible(s)
        if p:
            for n in p:
                stack.append((n, s[:], struc[:]))
        else:
            graphs.append(struc)

    return graphs


# --------------------------------------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    print('\nTesting connectivity')
    graph = RNAGraph(inlist=[0, 0, 1, 1, 2, 2])
    print('    ', graph.pairs)
    if not graph.connected():
        print('Not connected')

    print('\nlist format')
    structure = [0, 1, 2, 1, 0, 2]
    print('    input list', structure)
    graph = RNAGraph(inlist=structure)
    print('    serialized', str(graph))
    print('    pairs', graph.pairs)
    print('    list', graph.toList())
    print('    vienna', graph.toVienna())

    print('\n    reverse')
    print('    reversed', graph.reverse())
    print('    pairs', graph.pairs)
    print('    list', graph.toList())
    print('    vienna', graph.toVienna())

    print('\nenumerating: size, len, total')
    total = 0
    for size in range(1, 8):
        g = enumerateRNATopology(size)
        total += len(g)
        print('    ', size, len(g), total)

    print('\n3 stem graphs')
    graphs = enumerateRNATopology(3)
    for g in graphs:
        pgraph = RNAGraph(g)
        print('\ngraph:', g)
        print('    pairs:', pgraph.pairs, end='\t=>\t')
        print('    list:', pgraph.toList())
        print('    vienna', pgraph.toVienna())

        pgraph.reverse()
        print('    reversed:', str(pgraph))
        print('    pairs:', pgraph.pairs, end='\t=>\t')
        glist = pgraph.toList()
        print('    reversed list:', glist)
        print('    vienna', pgraph.toVienna())

    exit(0)
