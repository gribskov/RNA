class PairGraph:
    """=============================================================================================
    RNA graph class

    Synopsis
        graph = PairGraph()
        graph.fromList([0,1,1,0])
            or
        graph = PairGraph(l=[0,1,1,0])

        graph.reverse()     # reverse the order of stems left-right
    ============================================================================================="""

    def __init__(self, l=[]):
        """-----------------------------------------------------------------------------------------
        Graph constructor
        Internal data structure is a pair graph: a list in which each element is a stem.  The stems
        give the coordinates of the left and right half stems as integers along a line.
        -----------------------------------------------------------------------------------------"""
        self.pairs = []
        self.nstem = 0

        if l:
            self.fromList(l)

    def fromList(self, g):
        """"----------------------------------------------------------------------------------------
        convert a graph in array format to pair format.  returns a list of lists with the
        begin/end position of each stem
        :param g: graph in array format
        :return: number of stems in new graph
        -----------------------------------------------------------------------------------------"""
        self.nstem = int(len(g) / 2)
        self.pairs = [[] for i in range(self.nstem)]

        for i in range(len(g)):
            self.pairs[g[i]].append(i)

        return self.nstem

    def toList(self):
        """-----------------------------------------------------------------------------------------
        convert a graph in pair format to array format
        :return g: graph in list format
        -----------------------------------------------------------------------------------------"""
        g = [0 for i in range(self.nstem * 2)]

        stem = 0
        for pair in self.pairs:
            g[pair[0]] = stem
            g[pair[1]] = stem
            stem += 1

        return g

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

        return None


# --------------------------------------------------------------------------------------------------
# Non-object functions
# --------------------------------------------------------------------------------------------------
def possible(s):
    """
    return a list of the possible next stems.  Two kinds of stems can be added
    1) stems with use zero can be added in numerical order
    2) stems that have use == 1 can be added, but only if the the count of open stems would be > 0

    :param s: stem usage list
    :return: list of possible elements to add
    """
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


def enumerate(n):
    """---------------------------------------------------------------------------------------------
    enumerate all graphs with n stems
    :param n: number of stems
    :return: array of graph lists
    ---------------------------------------------------------------------------------------------"""
    graphs = []
    s = [0 for i in range(n)]
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

    graph = PairGraph(l=[0, 1, 2, 1, 0, 2])
    print(graph.pairs)
    graph.reverse()
    print(graph.pairs)

    total = 0
    for size in range(1, 8):
        g = enumerate(size)
        total += len(g)
        print(size, len(g), total)

    graphs = enumerate(3)
    for g in graphs:
        pgraph = PairGraph(l=g)
        print('\ngraph:', g)
        print('pairs:', pgraph.pairs, end='\t=>\t')
        print(pgraph.toList())
        pgraph.reverse()
        print('reversed:', pgraph.pairs)

    exit(0)
