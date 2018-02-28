class PairGraph:
    """=============================================================================================
    RNA graph class

    Synopsis
        from fasta import Fasta()

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


class Xios():
    """=============================================================================================
    Graph structure with XIOS edges
    ============================================================================================="""

    def __init__(self, pairs=''):
        """-----------------------------------------------------------------------------------------
        Xios constructor
        Initialize from graph
        :param graph: a PairGraph object
         -----------------------------------------------------------------------------------------"""
        self.edgelist = [[] for k in range(len(graph))]
        if pairs:
            self.getEdgeFromPairs(pairs)

        return None

    def __len__(self):
        """-----------------------------------------------------------------------------------------
        length is number of vertices

        :return: length of xios graph (number of vertices)
        -----------------------------------------------------------------------------------------"""
        return len(self.edgelist)

    def __str__(self):
        """-----------------------------------------------------------------------------------------
        formatted string representation of XIOS graph
        :return:
        -----------------------------------------------------------------------------------------"""
        s = '\n'
        for v in range(len(self.edgelist)):
            s += '{}:'.format(v)
            for e in self.edgelist[v]:
                s += ' {}{}'.format(e[0], e[1])
            s += '\n'

        return s

    def getEdgeFromList(self, graph):
        """-----------------------------------------------------------------------------------------
        Construct edgelist structure from pairGraph
        :param graph:
        :return: size of graph (number of stems)
        -----------------------------------------------------------------------------------------"""
        pairs = PairGraph(l=graph)
        self.getEdgeFromPairs(pairs)

        return len(self)

    def getEdgeFromPairs(self, pairs):
        """-----------------------------------------------------------------------------------------
        Construct edgelist structure from pairList. By the construction of the pair list, it is
        guaranteed that pairs[i][0] < pairs[j][0].

        :param list:
        :return: size of graph (number of stems)
        -----------------------------------------------------------------------------------------"""
        self.edgelist = [[] for k in range(len(pairs))]
        for i in range(len(pairs)):
            for j in range(i + 1, len(pairs)):
                if pairs.pairs[i][1] < pairs.pairs[j][0]:
                    self.edgelist[i].append([j, 's'])
                elif pairs.pairs[i][1] > pairs.pairs[j][1]:
                    self.edgelist[i].append([j, 'i'])
                else:
                    self.edgelist[i].append([j, 'o'])

        return len(self)

    def dfs(self):
        """-----------------------------------------------------------------------------------------
        Determine canonical graph labeling using gspan DFS algorithm
        :return:
        -----------------------------------------------------------------------------------------"""
        pass

        return


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
    print(str(graph))

    total = 0
    for size in range(1, 8):
        g = enumerate(size)
        total += len(g)
        print(size, len(g), total)

    graphs = enumerate(3)
    for g in graphs:
        pgraph = PairGraph(g)
        print('\ngraph:', g)
        print('pairs:', pgraph.pairs, end='\t=>\t')

        xios = Xios(pgraph)
        print('xios', str(xios))
        pgraph.reverse()
        print('reversed:', str(pgraph))
        print('pairs:', pgraph.pairs, end='\t=>\t')
        l = pgraph.toList()
        print('reversed list:', l)
        xios.getEdgeFromList(l)
        print('xios (from list):', str(xios))

    exit(0)
