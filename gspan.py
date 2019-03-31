import sys
import copy
import random
from functools import total_ordering


@total_ordering
class Edge(list):
    """=============================================================================================
    Edge class is necessary to implement lexicographic sorting
    for two edges a = (i1, j1, e1) and b = (i1, i2, e2), a < eb if
        i1 == i2 and j1 < j2        backward edge < forward edge
        i1 < j1  and j1 == i2       a is forward edge, b is an extension of j1

        edge ordering is transitive so if a < b and b < c, a < c

    ============================================================================================="""
    g2d = []  # class variable for translation of edges to dfs numbering

    def __init__(self, edge=None):
        """-----------------------------------------------------------------------------------------
        An edge is list of vertex 0, vertex 1, and edge type; all ints.  No assumption is made that
        v0 is less than v1.
        -----------------------------------------------------------------------------------------"""
        if edge is None:
            super(Edge, self).__init__([None, None, None])
        else:
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


# end of class Edge

class Gspan:
    """=============================================================================================
    Yan and Han algorithm for canonical graph labeling
    for now, assume a graph list of edges [v0, v1, e] where v0 and v1 are graph labels (arbitrary
    int) and e is the edge type (ordered: i, j, o, s )

    for a list of edges G
        init:
            sort G buy edgetype
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

    def graph_randomize(self):
        """-----------------------------------------------------------------------------------------
        randomly relaable the graph vertices.  this isw usful for generating graphs that have the
        same canonical form.

        :return: int, number of edges
        -----------------------------------------------------------------------------------------"""

        if self.graph is None:
            sys.stderr.write('Gspan.graph_randomize - graph is undefined\n')

        # make a dictionary of the vertex labels in the current graph
        vertex = {}
        v_original = []
        for edge in self.graph:
            for v in range(0, 2):
                if edge[v] in vertex:
                    vertex[edge[v]] += 1
                else:
                    vertex[edge[v]] = 1
                    v_original.append(edge[v])

        # new labels are sequential integers with randomized order, repeat until
        # an order different from the original is produced
        # adding extra makes sure the new vertices will not be contiguous beginning at zero
        extra = 4
        while True:
            v_new = [x for x in range(len(vertex) + extra)]
            random.shuffle(v_new)
            v_new = v_new[:-2]
            if not v_original == v_new:
                break

        # replace old labels with new labels and return the old -> new map as a dict
        i = 0
        for v in vertex:
            vertex[v] = v_new[i]
            i += 1

        for edge in self.graph:
            for v in range(0, 2):
                if edge[v] in vertex:
                    edge[v] = vertex[edge[v]]

            if edge[0] > edge[1]:
                # make i always less than j, flipping edgetype 0/1 if necessary
                edge[0], edge[1] = edge[1], edge[0]
                if edge[2] < 2:
                    edge[2] += 1
                    edge[2] = edge[2] % 2

        return vertex

    def save(self, edge, row):
        """-----------------------------------------------------------------------------------------
        push a partial solution onto the unexplored stack.
        Stores copy of g2d, edge, row

        :return: integer, length of stack
        -----------------------------------------------------------------------------------------"""
        g2d = copy.deepcopy(self.g2d)
        self.unexplored.append((g2d, edge, row))

        return len(self.unexplored)

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
        d_next = -1  # insures first edge is 0
        for d in self.g2d:
            if d is None:
                continue
            d_next = max(d, d_next)

        d_next += 1

        # update g2d, a restore is always an extension
        if self.g2d[edge[0]] is None:
            self.g2d[edge[0]] = d_next
            d_next += 1
        if self.g2d[edge[1]] is None:
            self.g2d[edge[1]] = d_next
            d_next += 1

        return d_next, row

    def sort(self, begin=0):
        """-----------------------------------------------------------------------------------------
        At each step, the new minimum edge in the unordered part of the graph is added to the DFS.
        This method sorts the unordered portion of the graph, beginning at row=begin, in DFS order
        Method:
            partition into three groups
                backward edges: both vertices known
                forward extension: one vertex  known
                both unknown
            copy the groups back into the graph sorted by appropriate vertices (varies by group)

        :return: graph
        -----------------------------------------------------------------------------------------"""
        graph = self.graph
        g2d = self.g2d

        # make separate lists of possible forward, backward, and unknown edges.  Makes sorting
        # easier

        backward = []
        forward = []
        unknown = []

        for edge in graph[begin:]:
            # the labels in the input graph, g, can be anything so they are treated as
            # the labels in the DFS code are integers or None

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

        # sort backward, forward, unknown and unknow edges and add to new order
        neworder = sorted(backward, key=lambda v: g2d[v[1]])

        forward.sort(key=lambda v: v[2])
        forward.sort(key=lambda v: g2d[v[0]], reverse=True)
        neworder += forward

        neworder += unknown

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
    graphset(0:2) have the same canonical representation as graphset[0]
    '''
    graphset = [[[0, 1, i], [1, 2, i], [2, 0, j]],
                [[0, 1, i], [0, 2, j], [1, 2, j]],
                [[0, 1, j], [0, 2, j], [1, 2, j]],
                [[0, 1, 1], [0, 2, 0], [0, 3, 0], [1, 2, 0], [1, 3, 0], [2, 3, 2]],
                [[0, 1, 0], [0, 2, 0], [0, 3, 0], [0, 4, 0], [1, 4, 2], [2, 3, 0], [2, 4, 2],
                 [3, 4, 2]]]

    # graph normalization create an unnormalized graph by doubling the vertex numbers

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

    for graph in graphset:
        print('\nGraph normalization')
        g = copy.deepcopy(graph)
        print('    original graph: {}'.format(g))

        # for edge in g:
        #     for i in range(0, 2):
        #         edge[i] *= 2
        # print('    un-normalized graph: {}'.format(g))

        gspan = Gspan(graph=g)
        map = gspan.graph_randomize()
        # print('    map', map)
        print('    randomized graph: {}'.format(gspan.graph))
        gspan.graph_normalize()
        print('    renormalized graph: {}'.format(gspan.graph))

    # g = graphset[1]
    # g = [[0, 1, 1], [0, 2, 0], [0, 3, 0], [1, 2, 0], [1, 3, 0], [2, 3, 2]]
    print('\nGspan canonical graph')
    g = graphset[4]
    print('\n\tinput graph', g)
    gspan = Gspan(graph=g)
    # map = gspan.graph_randomize()
    gspan.graph_normalize()
    print('\trenormalized graph: {}'.format(gspan.graph))
    glen = len(gspan.graph)

    # initialize first edge after sorting by edgetype = e[2]: i < j < o < s < x
    gspan.graph = sorted(gspan.graph, key=lambda e: e[2])

    row = 0
    e_first = gspan.graph[0][2]
    for edge in gspan.graph:
        if edge[2] == e_first:
            gspan.save(edge, row)

    # TODO add to min DFS

    while gspan.unexplored:
        d, row = gspan.restore()
        print('\n\trestored d={} row={} edge={}'.format(d, row, gspan.graph[row]))
        if d is None:
            break

        g2d = gspan.g2d
        row += 1

        while row < glen:
            # sort the possible extension by backward, forward, unknown, the current row is the
            # best extension (sorted first in list
            gspan.sort(begin=row)
            edge = gspan.graph[row]

            print('\n\tb graph', gspan.graph, '\n\t\tdfs', gspan.graph2dfs(), '\n\t\tg2d',
                  gspan.g2d)

            while row < glen and g2d[edge[1]] is not None:
                # add all backward edges, they are always unique and never require resorting
                row += 1
                if row >= glen:
                    break
                edge = gspan.graph[row]

            if row < glen:
                # forward extension, there must be one or we are at the end of the graph
                # gspan.g2d[edge[1]] = d
                e_first = edge[2]
                v0_first = gspan.graph[row][0]

                # if there are equivalent extensions, save them
                # v0 defined, v1 undefined, edgetype = e_first edgetype
                rr = row + 1
                while rr < glen \
                        and gspan.graph[rr][0] == v0_first \
                        and gspan.graph[rr][2] == e_first:
                    gspan.save(gspan.graph[rr], row)
                    rr += 1

                gspan.g2d[edge[1]] = d
                d += 1
                row += 1

        # end of loop over rows of dfs code

        # print('\ngraph', gspan.graph, '\n    dfs', gspan.graph2dfs(), '\n    g2d', gspan.g2d)
        print('\t----------------------------------------')

    # end of loop over all starting vertices

exit(0)
