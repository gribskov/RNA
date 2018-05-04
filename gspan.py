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
        return (self[0] == other[0] and self[1] == other[1] and self[2] == other[2])

    def __ne__(self, other):
        return not (self[0] == other)

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

    def flip(self, row=0):
        """-----------------------------------------------------------------------------------------
        convert all j edges to i edges from row to end of graph
        :param row: integer, beginning row
        :return: integer, number flipped
        -----------------------------------------------------------------------------------------"""
        n = 0
        for i in range(row, len(self.graph)):
            if self.graph[i][2] == 1:
                self.graph[i].reverse()
                n += 1

        return n

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

    def restore(self):
        """-----------------------------------------------------------------------------------------
        pop a saved dfs code from the unexplored stack
        swap the edge with the current edge at the specified row
        flip j edges in rows following the new edge
        :return: integer, next available dfs vertex
        -----------------------------------------------------------------------------------------"""
        if len(self.unexplored) == 0:
            return None, None

        g2d, edge, row = self.unexplored.pop()
        self.g2d = g2d

        epos = self.graph.index(edge)
        self.graph[row], self.graph[epos] = self.graph[epos], self.graph[row]
        self.flip(row)

        all_undef = True
        for d in self.g2d:
            if d is None:
                continue
            d_next = max(d, d_next)
            all_undef = False
        if all_undef:
            d_next = 0
        else:
            d_next += 1

        # update g2d, restore are always an extension
        if self.g2d[edge[0]] is None:
            self.g2d[edge[0]] = d_next
            d_next += 1
        if self.g2d[edge[1]] is None:
            self.g2d[edge[1]] = d_next
            d_next += 1

        return d_next, row

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
                    # both i and j are none, do nothing
                    continue
                else:
                    # i is none, j is known, reverse the edge
                    edge.reverse()

            elif self.g2d[edge[1]] is None:
                # edge 0 is not None, tested above
                continue

            elif self.g2d[edge[0]] < self.g2d[edge[1]]:
                # edge.reverse()
                continue

        return True

    def set_i_edge(self):
        """-----------------------------------------------------------------------------------------
        When a new first edge is popped from the stack, flip all j edges back to i edges
        :return: True
        -----------------------------------------------------------------------------------------"""
        for edge in self.graph:
            if edge[2] == 1:
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

    def sort(self, begin=0):
        """-----------------------------------------------------------------------------------------
        At each step, the new minimumedge in the unordered part of the graph is added to the DFS.
        This method sorts the unordered portion of the graph, beginning at row=begin, in DFS order
        Method:
            partition into three groups
                backward edges: both vertices known
                forward extension: one vertex  known
                both unknown
            copy the groups back into the graph sorted by appropriate vertices (varies by group)

        :return:
        -----------------------------------------------------------------------------------------"""
        graph = self.graph
        g2d = self.g2d

        backward = []
        forward = []
        unknown = []
        for edge in graph[begin:]:
            if g2d[edge[0]] is None:
                # vertex 0 undefined
                if g2d[edge[1]] is None:
                    # both undefined: can be sorted by edgetype only
                    # should not need to flip, normalize does it
                    unknown.append(edge)
                else:
                    # vertex 0 undefined, vertex 1 defined: possible extension
                    # defined vertex should be v0, so reverse the edge
                    edge.reverse()
                    forward.append(edge)
            else:
                # vertex 0 defined
                if g2d[edge[1]] is None:
                    # vertex 0 defined, vertex 1 undefined: possible forward extension
                    forward.append(edge)
                else:
                    # both defined, backward edge
                    # v0 must be > v1, so reverse edge if needed
                    if g2d[edge[0]] < g2d[edge[1]]:
                        edge.reverse()
                    backward.append(edge)

        # copy edges into graph: backward, forward, unknown
        neworder = []
        for edge in sorted(backward, key=lambda v: g2d[v[1]]):
            # backward edges should only come from the rightmost vertex and are sorted by v1
            neworder.append(edge)

        for edge in sorted(forward, key=lambda v: g2d[v[0]], reverse=True):
            # forward extensions are made from the largest v0 first
            neworder.append(edge)

        for edge in unknown:
            neworder.append(edge)

        graph[begin:] = neworder

        return graph

    # end of sort

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
    # these tests don't make sense
    #  if g == gspan.graph:
    #     print('      passes test')
    # else:
    #     print('      fails test')
    # gspan.g2d = [1, 0, 2, 3]
    # gspan.sort_by_edge()
    # print('    sorted graph using g2d={}: {}\t{}'.format(gspan.g2d, gspan.graph,gspan.graph2dfs()))

    print('\nEdge manipulation\n')
    e = Edge()
    e.set(2, 3, 0)
    e.g2d = [0, 1, 2]
    print('    edge', e)
    e.reverse()
    print('    edge reversed', e)
    e.set(1, 2, 1)
    e.g2d = [2, 1, 0]
    print('    dfs numbering using {}: {}'.format(e.g2d, e))
    e.reverse()
    print('    dfs numbering reversed', e)

    print('Edge comparison\n')
    Edge.g2d = [1, 0, 2, None]
    e1 = Edge([0, 3, 0])
    e2 = Edge([2, 1, 1])
    print('    edge 1', e1)
    print('    edge 2', e2)

    if e1 < e2:
        print('e1 smaller')
    if e2 < e1:
        print('e2 smaller')
    g = [Edge([0, 1, 0]), Edge([2, 0, 1])]
    Edge.g2d = [0, None, 2, 3]

    # g = graphset[1]
    g = [[0, 1, 1], [0, 2, 0], [0, 3, 0], [1, 2, 0], [1, 3, 0], [2, 3, 2]]
    print('\ninput graph', g)
    gspan = Gspan(graph=g)
    glen = len(gspan.graph)

    # initialize first edge
    # gspan.sort_by_edge()
    graph = gspan.graph

    row = 0
    first_e = graph[0][2]
    for edge in graph:
        if edge[2] == first_e:
            g2d = copy.deepcopy(gspan.g2d)
            gspan.unexplored.append((g2d, edge, row))

    while gspan.unexplored:
        d, row = gspan.restore()
        g2d = gspan.g2d
        if d is None:
            break

        row += 1

        while row < glen:
            gspan.sort(begin=row)
            print('\nb graph', gspan.graph, '\n    dfs', gspan.graph2dfs(), '\n    g2d', gspan.g2d)
            edge = gspan.graph[row]
            # add all backward edges, they are always unique and never require sorting
            while row < glen and g2d[edge[1]] is not None:
                row += 1
                if row >= glen:
                    break
                edge = gspan.graph[row]

            if row < glen and g2d[edge[1]] is None:
                # forward extension
                gspan.g2d[edge[1]] = d
                d += 1
                row += 1

        # end of loop over rows of dfs code

        print('\ngraph', gspan.graph, '\n    dfs', gspan.graph2dfs(), '\n    g2d', gspan.g2d)
        print('----------------------------------------')

        # exit(0)

    # end of loop over all starting vertices

exit(0)
