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
        return self[0] == other[0] and self[1] == other[1] and self[2] == other[2]

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

    The input graph G can have any arbitrary labels (strings).
        First it must be converted to a normalized graph where the vertices are sequential integers
        for testing, a graph may be entered as a python array of arrays (see from_list)

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
        self.graph = []         # graph in normalized labelling
        self.map = None
        self.vnum = 0           # number of vertices in graph
        self.vnext = 0          # next dfs vertex available to use
        self.mindfs = []
        self.g2G = []           # list to convert g indices to original labels
        self.g2d = []           # list to covert g labels do dfs labels
        self.d2g = []           # list to conver dfs labels to g labels
        self.unexplored = []    # stack of partial solutions that need to be searched
                                # [d2g, edge, row_num]

        if graph:
            if type(graph) is list:
                self.from_list(graph)
                self.graph_normalize()

            elif type(graph) is str:
                self.from_string(graph)

            else:
                sys.stderr.write('Gspan::__init__ - unknown graph type ({})\n'.format(graph))

    def from_list(self, graph):
        """-----------------------------------------------------------------------------------------
        load a graph in the form of a list, into a graph object composed of Edge()

        :param graph: list of lists:
        :return: int, number of vertices
        -----------------------------------------------------------------------------------------"""
        self.graph = list()
        v = 0
        for edge_in in graph:
            edge = Edge(edge_in)
            self.graph.append(edge)
            v += 1

        self.mindfs = [Edge() for _ in self.graph]

        self.vnum = v
        return v

    def from_string(self, graphstr):
        """-----------------------------------------------------------------------------------------
        Reads a graph from a string and converts to non_normalized integers.  For now graphs are
        a list of edges, each edge is a triple of vertex0 vertex1 edge_type.  All are abritrary
        strings without quotation marks. square brackets and commas are ignored so a python list of
        lists is OK

        :param graphstr:
        :return:
        -----------------------------------------------------------------------------------------"""
        translation = str.maketrans('[],', '   ')
        graphstr = graphstr.translate(translation)
        elist = graphstr.split()
        G2g = {}  # hash showing translation of original labels to ints
        g2G = []  # back translate from g index to original labels
        g = []  # transformed graph
        v = 0
        for i in range(0, len(elist), 3):
            if elist[i] not in G2g:
                G2g[elist[i]] = v
                g2G.append(elist[i])
                v += 1

            if elist[i + 1] not in G2g:
                G2g[elist[i + 1]] = v
                g2G.append(elist[i + 1])
                v += 1

            v0 = G2g[elist[i]]
            v1 = G2g[elist[i + 1]]
            etype = int(elist[i + 2])
            self.graph.append([v0, v1, etype])

        self.g2G = g2G
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
        self.unexplored.append([g2d, edge, row])

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

        # all_undef = True
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

        unknown.sort(key=lambda v: v[2])

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

    def add_forward(self, row):
        """-----------------------------------------------------------------------------------------
        Add forward edge. should work also for first edge
            update g2d
            update d2g
            update highest used v

        :param row:
        :return:
        -----------------------------------------------------------------------------------------"""
        # TODO need to implement
        pass

        return True

    def isbackward(self, row):
        """-----------------------------------------------------------------------------------------
        True if the edge specified by row is backward. for backward edges, v0 > v1 in dfs labeling
        for a forward edge, only one vertex is defined. if both are defined, it should be a backward
        edge
        TODO: if this is true the test can be simpler


        :param row: int, row of the dfs code
        :return: logical
        -----------------------------------------------------------------------------------------"""
        # edge = self.graph[row]
        g2d0= self.g2d[self.graph[row][0]]
        g2d1 = self.g2d[self.graph[row][1]]
        if g2d0 and g2d1:
            if g2d0 > g2d1:
                return True

        return False


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

    # testing reading graphs from string
    gstr = '[[0, 1, 1], [0, 2, 0], [0, 3, 0], [1, 2, 0], [1, 3, 0], [2, 3, 2]]'
    g = Gspan(gstr)
    gstr = '[[a, c, 1], [a, b, 0], [a, d, 0], [c, b, 0], [c, d, 0], [b, d, 2]]'
    g = Gspan(gstr)

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

    # beginning of Gspan algorithm

    gspan.sort()
    row = 0

    # add all edges that are the sames type as the first edge to the stack
    first_edge_type = gspan.graph[0][2]
    for edge in gspan.graph:
        if edge[2] == first_edge_type:
            gspan.save(edge, row)

    while gspan.unexplored:
        d, row = gspan.restore()
        print('\n\trestored d={} row={} edge={}'.format(d, row, gspan.graph[row]))
        if d is None:
            # not sure when this is supposed to happen?
            break

        g2d = gspan.g2d
        row += 1

        # sort the possible extension by backward, forward, unknown, the current row is the
        # best extension (sorted first in list)
        # TODO check that backward edges don't need to be added when restoring
        # i think they do
        gspan.sort(begin=row)

        while row < glen:

            # the first edge should always be a forward extension because we add any backward
            # extensions at the bottom of this loop.  backwards edges never need resorting
            # because no new vertices are added
            # The next edge should never be an unknown edge because that would indicate two
            # disjoint graphs are present
            edge = gspan.graph[row]

            print('\n\tb graph', gspan.graph, '\n\t\tdfs', gspan.graph2dfs(), '\n\t\tg2d',
                  gspan.g2d)

            # forward extension:
            # the next sorted edge is a forward extension, d2g and g2d need update
            # any equivalent forward extensions (same v0, seme edge tpe, have to go on stack)
            edge_type = edge[2]
            v0_first = gspan.graph[row][0]
            gspan.add_forward(row)

            rr = row + 1
            while rr < glen \
                    and gspan.graph[rr][0] == v0_first \
                    and gspan.graph[rr][2] == edge_type:
                gspan.save(gspan.graph[rr], row)
                rr += 1

            # add all backward edges, no new sorting is required

            row += 1
            while gspan.isbackward(row):
                edge = gspan.graph[row]
                row += 1
                if row >= len(gspan.graph):
                    break
            # check dfs

            # should now be at the non-backward edge, there are three possibilities
            #   another forward edge on rightmost path
            #   and unknown edge (needs resort)
            #   or the end of the list

            if row == glen:
                # graph is done
                break

            gspan.sort(begin=row)

        # end of loop over rows of dfs code

        # print('\ngraph', gspan.graph, '\n    dfs', gspan.graph2dfs(), '\n    g2d', gspan.g2d)
        print('\t----------------------------------------')

    # end of loop over all starting vertices

exit(0)
