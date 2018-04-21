import sys
import copy
from functools import total_ordering


@total_ordering
class Edge(list):
    """=============================================================================================
    Edge class is necessary to implement lexicographic sorting
    e1 < e2 if i1==i2 and j1<j2
               i1<j1 and j1=i2
               e1 <(E,T) and e2 <(E,T) e3 then e1 < e3
    ============================================================================================="""
    g2d = []  # class variable for translation of edges to dfs numbering

    def __init__(self, edge=[None, None, None]):
        """-----------------------------------------------------------------------------------------
        An edge is list of vertex 0, vertex 1, and edge type; all ints.  No assumption is made that
        v0 is less than v1.
        -----------------------------------------------------------------------------------------"""
        super(Edge, self).__init__(edge)

    def __eq__(self, other):
        return ( self[0] == other[0] and self[1] == other[1] and self[2] == other[2] )

    def __ne__(self, other):
        return not ( self[0] == other )

    def __lt__(self, other):
        ia = Edge.g2d[self[0]]
        ja = Edge.g2d[self[1]]

        ib = Edge.g2d[other[0]]
        jb = Edge.g2d[other[1]]

        if ia is None:
            # a is unmapped
            if ib is None:
                # b is unmapped, smaller edgetype is less
                return self[2] < other[2]
            else:
                # a unmapped, b is known or on rightmost path
                return False
        elif ib is None:
            # a is known or on rightmost path, b is unmapped
            return True

        # both a and b are at least partially known
        if ja is None:
            # a is an extension
            if jb is None:
                # a and b are both extensions on rightmost path
                return ia > ib
            else:
                # a is extension, b is known
                return False
        elif jb is None:
            # a is known, b is extension
            return True

        # both a and b are known edges
        if ia < ja:
            # a is forward
            if ib > jb:
                # b is backward
                return False
            else:
                # both forward edges
                if ja < jb:
                    return True
                if ja > jb:
                    return False

                # ja = jb
                if ia < ib:
                    return True
                if ia > ib:
                    return False

        else:
            # a is backward
            if ib > jb:
                # both backward edges
                if ia < ib:
                    return True
                if ia > ib:
                    return False

                # ia ==ib
                if ja < jb:
                    return True
                return False

            else:
                # a backward, b forward
                return True

    def set(self, v0=None, v1=None, e=None):
        """-----------------------------------------------------------------------------------------
        set v0, v1, edge
        :return: True
        -----------------------------------------------------------------------------------------"""
        self[0] = v0
        self[1] = v1
        self[2] = e

        return True

    def reverse(self):
        """-----------------------------------------------------------------------------------------
        reverse the direction of the edge

        :return: True
        -----------------------------------------------------------------------------------------"""
        self[0], self[1] = self[1], self[0]
        if self[2] < 2:
            self[2] ^= 1

        return True


class Gspan:
    """=============================================================================================
    Yan and Han algorithm for canonical graph labeling
    for now, assume a graph list of edges [v0, v1, e] where v0 and v1 are graph labels (arbitrary
    int) and e is the edge type (ordered: i, j, o, s )

    for a list of edges G
        init:
            sort G
            push G for all equal edges on unexplored
            
        main loop:
            pop G from unexplored
            add to dfs
            test vs mindfs
            
            update g2d (graph to dfs labelling) map
            sort G
            add all backward edges to dfs
            push equivalent forward extensions on unexplored
            
    stack has G, d2g and ?

    much better
    all you need to save is g2d list and row of edgelist
    sorted edgelist is the dfs
            
    Michael Gribskov     20 April 2018
   
    ============================================================================================="""

    def __init__(self, graph=None):
        """-----------------------------------------------------------------------------------------
        gspan constructor
        -----------------------------------------------------------------------------------------"""
        self.graph = None
        self.map = None
        self.vnum = 0
        self.mindfs = []
        self.g2d = None
        self.d2g = None
        self.unexplored = []

        if graph:
            self.graph_load(graph)
            self.graph_normalize()

    def graph_load(self, graph):
        """-----------------------------------------------------------------------------------------
        load a graph in the form of a list, into a graph object composed of Edge()

        :param graph: list of lists
        :return: int, number of vertices
        -----------------------------------------------------------------------------------------"""
        self.graph = list()
        v = 0
        for edge_in in graph:
            edge = Edge(edge_in)
            self.graph.append(edge)
            v += 1

        self.vnum = v
        return v

    def graph_normalize(self):
        """-----------------------------------------------------------------------------------------
        Count the number of vertices and, if necessary, renumber so they begin at zero
        initialize d2g and g2d to be the correct size
        :return: int, number of vertices
        -----------------------------------------------------------------------------------------"""
        if self.graph is None:
            sys.stderr.write('Gspan.graph_normalize - graph is undefined\n')

        # list of vertices, index in v gives the normalized index
        v = []
        for edge in self.graph:
            for i in range(0, 2):
                if edge[i] not in v:
                    v.append(edge[i])

        # convert edge numbers, flip directed edges so they are all 0 (i) not 1 (j)
        for edge in self.graph:
            for i in range(0, 2):
                edge[i] = v.index(edge[i])
            if edge[2] == 1:
                edge.reverse()

        self.vnum = len(v)
        self.map = v

        # initialize d2g and g2d
        self.d2g = [None for _ in range(0, self.vnum)]
        self.g2d = [None for _ in range(0, self.vnum)]

        return self.vnum

    def order(self, row=0):
        """-----------------------------------------------------------------------------------------
        given the g2d map, orient the edges so that when only one vertex is defined it is v0 (i)
        or when unadded known edges have v0 < v1
        :param row: int, beginning row for transformation
        :return: True
        -----------------------------------------------------------------------------------------"""
        for edge in self.graph[row:]:
            if self.g2d[edge[0]] is None:
                if self.g2d[edge[1]] is None:
                    continue
                edge.reverse()

            elif self.g2d[edge[1]] is None:
                    continue

            elif self.g2d[edge[0]] < self.g2d[edge[1]]:
                edge.reverse()

        return True

    def sort_by_edge(self, row=0):
        """-----------------------------------------------------------------------------------------
        sort graph by edges, using current g2d map
        :return:
        -----------------------------------------------------------------------------------------"""
        Edge.g2d = self.g2d
        tmp = self.graph[row:]
        tmp.sort()
        self.graph[row:] = tmp
        print('    sorted:', self.graph)

        return True

    def graph2dfs(self):
        """-----------------------------------------------------------------------------------------
        convert the node labels in the graph to dfs labels
        :return: list of edges
        -----------------------------------------------------------------------------------------"""
        g2d = self.g2d
        dfs = []
        for edge in self.graph:
            row = []
            dfs.append(row)
            for i in range(0, 2):
                row.append(g2d[edge[i]])
            row.append(edge[2])

        return dfs


# ==================================================================================================
# testing
# ==================================================================================================
if __name__ == '__main__':
    # mapping of edge types to integers to avoid quoting all the time
    i = 0
    j = 1
    o = 2
    s = 3

    '''
    graphset(0:2) have the canonical representation graphset[0]
    '''
    graphset = [[[0, 1, i], [1, 2, i], [2, 0, j]],
                [[0, 1, i], [0, 2, j], [1, 2, j]],
                [[0, 1, j], [0, 2, j], [1, 2, j]],
                [[0, 1, 1], [0, 2, 0], [0, 3, 0], [1, 2, 0], [1, 3, 0], [2, 3, 2]]]

    # graph normalization create an unnormalized graph by doubling the vertex numbers
    print('Graph normalization')
    g = copy.deepcopy(graphset[1])
    print('    original graph: {}'.format(g))
    for edge in g:
        for i in range(0, 2):
            edge[i] *= 2
    print('    un-normalized graph: {}'.format(g))

    gspan = Gspan(graph=g)
    # graph normalization should be automatic
    # gspan.graph_normalize()
    print('    renormalized graph: {}'.format(gspan.graph))
    if g == gspan.graph:
        print('    passes test')

    e = Edge()
    e.set(2, 3, 0)
    e.g2d = [0, 1, 2]
    print(e)
    e.reverse()
    print(e)
    e.set(1, 2, 1)
    e.g2d = [2, 1, 0]
    print(e)
    e.reverse()
    print(e)

    Edge.g2d = [1, 0, 2, None]
    e1 = Edge([0, 3, 0])
    e2 = Edge([2, 1, 1])
    if e1 < e2:
        print('e1 smaller')
    if e2 < e1:
        print('e2 smaller')
    g=[ Edge([0,1,0]), Edge([2,0,1])]
    Edge.g2d=[0, None, 2, 3]
    g.sort()
    print(g)

    # g = graphset[1]
    g = [[0, 1, 1], [0, 2, 0], [0, 3, 0], [1, 2, 0], [1, 3, 0], [2, 3, 2]]
    print('input graph', g)
    gspan = Gspan(graph=g)
    glen = len(gspan.graph)

    # initialize first edge
    gspan.sort_by_edge()
    e = gspan.graph[0]
    gspan.g2d[e[0]] = 0
    gspan.g2d[e[1]] = 1
    d = 2
    i = 1

    while i < glen - 1:
        print('\ngraph', gspan.graph, '\n    dfs', gspan.graph2dfs(), '\n    g2d', gspan.g2d)
        gspan.order(i)
        print('    order', gspan.graph)
        gspan.sort_by_edge(i)
        print('    sorted dfs', gspan.graph2dfs())
        e = gspan.graph[i]
        if gspan.g2d[e[1]] is None:
            # forward extension
            gspan.g2d[e[1]] = d
            d += 1

        i += 1

    print('\ngraph', gspan.graph, '\n    dfs', gspan.graph2dfs(), '\n    g2d', gspan.g2d)

exit(0)
