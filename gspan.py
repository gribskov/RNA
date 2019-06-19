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

    def copy(self):
        return Edge(self[:])

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


    partial solution stack: gspan.unexplored stores
    g2d list, current edge, and row number

    Michael Gribskov     20 April 2018
    ============================================================================================="""

    def __init__(self, graph=None):
        """-----------------------------------------------------------------------------------------
        gspan constructor
        -----------------------------------------------------------------------------------------"""
        self.graph = []  # graph in normalized labelling
        self.nforward = 0  # number of forward edges in sorted graph
        self.nbackward = 0  # number of backward edges in sorted graph
        self.nunknown = 0  # number of unknown edges in sorted graph
        # self.map = None
        self.vnum = 0  # number of vertices in graph
        self.vnext = 0  # next dfs vertex available to use
        self.mindfs = []
        self.mindfslen = 0
        self.g2G = []  # list to convert g laberls (indices) to original labels
        self.g2d = []  # list to covert g labels to dfs labels
        self.d2g = []  # list to conver dfs labels to g labels
        self.row = 0  # row in current dfs
        self.unexplored = []  # stack of partial solutions that need to be searched
        # [d2g, edge, row_num]

        if graph:
            if isinstance(graph, list):
                self.from_list(graph)

            elif isinstance(graph, str):
                self.from_string(graph)

            else:
                sys.stderr.write('Gspan::__init__ - unknown graph type ({})\n'.format(graph))

            self.graph_normalize()

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
        # g = []  # transformed graph
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
        TODO check if this is the same as Edge.reverse()

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
        #        self.map = v

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
        push a partial solution onto the unexplored stack. Stores copy of g2d, edge, row

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
            return False

        g2d, edge, row = self.unexplored.pop()
        self.g2d = g2d

        if row == 0:
            if edge[2] == 1:
                edge.reverse()
        elif g2d[edge[0]] is None:
            # next edge is always a forward extension. if v0 is not defined in g2d, flip edge
            edge[0], edge[1] = edge[1], edge[0]
            if edge[2] < 2:
                edge[2] ^= 1

        # this makes the popped edge the next to be added
        try:
            epos = self.graph.index(edge)
        except ValueError:
            epos = self.graph.index([edge[1],edge[0],edge[2]])
        self.graph[epos] = self.graph[row]
        self.graph[row] = edge
        # self.flip(row)

        # rebuild d2g from g2d
        self.d2g = [None for _ in self.d2g]
        self.vnext = 0
        for i in range(len(self.g2d)):
            if self.g2d[i] is not None:
                self.d2g[self.vnext] = i
                self.vnext += 1

        # add the saved edge
        # if self.g2d[edge[0]] is None:
        #     # edge is backward, v0 always defined for an extension
        #     edge.reverse()

        if self.vnext == 0:
            self.g2d[edge[0]] = 0
            self.d2g[0] = edge[0]
            self.vnext = 1

        self.g2d[edge[1]] = self.vnext
        self.d2g[self.vnext] = edge[1]
        self.vnext += 1

        self.row = row + 1
        return True

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

        self.nforward = len(forward)
        self.nbackward = len(backward)
        self.nunknown = len(unknown)

        # sort backward, forward, unknown and unknown edges and add to new order

        # backward edges are sorted by v1.  v0 should all be the same
        neworder = sorted(backward, key=lambda v: g2d[v[1]])

        # forward edges have defined v0, and undefined v1.  sort by edge type (v2)
        # the largest v0 is the rightmost edge, smaller v0 are extensions on the rightmost path
        forward.sort(key=lambda v: v[2])
        forward.sort(key=lambda v: g2d[v[0]], reverse=True)
        neworder += forward

        # edges with neither vertex defined are sorted by edge type(v2)
        unknown.sort(key=lambda v: v[2])

        neworder += unknown

        graph[begin:] = neworder

        return graph

    # end of sort

    def graph2dfs(self, last):
        """-----------------------------------------------------------------------------------------
        convert the node labels in the graph to dfs labels

        :last: integer, last row to copy
        :return: list of edges
        -----------------------------------------------------------------------------------------"""
        g2d = self.g2d
        dfs = []
        for i in range(last):
            self.mindfs[i] = [g2d[self.graph[i][0]], g2d[self.graph[i][1]], self.graph[i][2]]

        self.mindfslen = last

        return dfs

    def minDFS(self):
        """-----------------------------------------------------------------------------------------
        reimplement gspan

        :return:
        -----------------------------------------------------------------------------------------"""
        mindfs = []
        row = 0
        searching = True

        while searching:
            # sort
            self.sort(begin=row)

            # add backward edges, since the graph is sorted, all this requires is updating the
            # begin point in the graph array
            row += self.nbackward

            if row < len(self.graph):
                # skip adding forward edges if done

                if row == 0 and self.graph[0][2] == 2:
                    # special case when all edges are undirected
                    edge = self.graph[0].copy()
                    edge.reverse()
                    self.save(edge, row)
                # save equivalent forward edges.  We know the number of forward edges from the sort,
                # the edges are equivalent if they have the same v0 and edge type (v2)
                # when there are zero forward edges, which should only occur for the first edge,
                # the equivalent unknown edges should be added since v0 = None and edge=minedge
                first_edge_type = self.graph[row][2]
                v0 = self.g2d[self.graph[row][0]]
                for edge in self.graph[row + 1:]:
                    if edge[2] != first_edge_type:
                        break
                    if self.g2d[edge[0]] != v0:
                        break

                    # if you pass these tests, fall though to saving this edge on stack
                    self.save(edge, row)
                    if row == 0 and edge[2] == 2:
                        # special case when all edges are undirected
                        erev = edge.copy()
                        erev.reverse()
                        self.save(erev, row)

                if row == 0:
                    # add the first edge
                    self.d2g[0] = self.graph[0][0]
                    self.d2g[1] = self.graph[0][1]
                    self.vnext = 2
                    self.g2d[self.d2g[0]] = 0
                    self.g2d[self.d2g[1]] = 1
                else:
                    self.d2g[self.vnext] = self.graph[row][1]
                    self.g2d[self.graph[row][1]] = self.vnext
                    self.vnext += 1

                row += 1

            # check for minimum dfs
            if not self.minimum(row) or row == len(self.graph):
                # print('current minDFS', self.mindfs)
                # print('         graph', self.graph, '\n')
                searching = self.restore()
                row = self.row

        return self.mindfs

    def minimum(self, row):
        """-----------------------------------------------------------------------------------------
        check if the current dfs is <= the min dfs.
        mindfs is stored in d space, current dfs must be compared in d space

        :return: logical
        -----------------------------------------------------------------------------------------"""

        if row > self.mindfslen:
            self.graph2dfs(row)
            return True

        cmp = None
        for i in range(row):
            cedge = self.edge_g2d(self.graph[i])
            medge = self.mindfs[i]

            cdir = self.edge_dir(cedge)
            mdir = self.edge_dir(medge)

            cmp = 'eq'
            if cdir < mdir:
                # c backward, m forward, current is lt
                cmp = 'lt'
                break

            if cdir > mdir:
                # c forward, m backward, min is lt
                cmp = 'gt'
                break

            # directions are equal
            if cdir is 'f':
                # forward edges
                if cedge[0] > medge[0]:
                    # c is lt
                    cmp = 'lt'
                    break
                elif cedge[0] < medge[0]:
                    # m is lt
                    cmp = 'gt'
                    break
                else:
                    # equal v0, v1 must always be the same
                    if cedge[1] != medge[1]:
                        sys.stdout.write('gspan::minimum - forward v1 not equal\n')
                    if cedge[2] < medge[2]:
                        # c is lt
                        cmp = 'lt'
                        break
                    elif cedge[2] > medge[2]:
                        # m is lt
                        cmp = 'gt'
                        break
                    else:
                        # edges are equal
                        continue

            if cdir is 'b':
                # backward edges, v0 must be the same
                if cedge[0] != medge[0]:
                    sys.stdout.write('gspan::minimum - backward v0 not equal\n')

                if cedge[1] < medge[1]:
                    # c is lt
                    cmp = 'lt'
                    break
                elif cedge[0] > medge[0]:
                    # m is lt
                    cmp = 'gt'
                    break
                else:
                    # equal v0, check edge type
                    if cedge[2] < medge[2]:
                        # c is lt
                        cmp = 'lt'
                        break
                    elif cedge[2] > medge[2]:
                        # m is lt
                        cmp = 'gt'
                        break
                    else:
                        # edges are equal
                        continue

            # end of loop over edges

        if cmp is 'lt':
            # current is definitively less than minimum, save as new minimum
            self.graph2dfs(row)
            return True

        elif cmp is 'eq':
            return True

        # gt, fall through
        # print('not minimal\ncurrent minDFS', self.mindfs[:row - 1],
        #       self.edge_g2d(self.graph[row - 1]))
        # print('         graph', self.graph[:row], '\n')
        return False

    def edge_dir(self, edge):
        """-----------------------------------------------------------------------------------------
        determine whether edge is forward or backward. edge should be in d space, use edge_g2d() to
        convert if necessary.

        :param edge: Edge
        :return: string, 'f' or 'b'
        -----------------------------------------------------------------------------------------"""
        direction = 'b'
        if edge[0] < edge[1]:
            direction = 'f'

        return direction

    def edge_g2d(self, edge):
        """-----------------------------------------------------------------------------------------
        convert edge from g space to d space
        
        :param edge: Edge (g space)
        :return: Edge (d space)
        -----------------------------------------------------------------------------------------"""
        g2d = self.g2d
        return Edge([g2d[edge[0]], g2d[edge[1]], edge[2]])


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
    graphset = [
        # [[0,1,0], [0,2,0], [1,2,0], [0,3,2], [1,3,2], [2,3,0], [3,4,0]],
        # [[0,1,0], [0,2,2], [0,4,2], [1,2,2], [2,3,2], [2,4,0], [2,5,2], [4,5,2], [5,3,0]],
        # [[0,1,0], [0,2,0], [0,3,0], [0,4,0], [0,5,0], [1,2,0], [4,5,0]],
        # [[0, 1, 0], [0, 2, 0], [0, 3, 0], [0, 4, 0], [0, 5, 0], [1, 2, 0], [4, 5, 0], [0,6,0], [1,6,2], [4,6,2]],
        # [[0, 1, 0], [0, 2, 0], [0, 3, 0], [0, 4, 0], [0, 5, 0], [1, 2, 0], [4, 5, 0], [0, 6, 0],
        # [1, 6, 2], [4, 6, 2], [0,7,0], [1,7,2], [2,7,2], [4,7,2], [5,7,2]],
        # [[0, 1, 0], [0, 2, 0], [0, 3, 0], [0, 4, 0], [1, 4, 0], [0, 5, 0], [2, 5, 0], [4, 6, 2],
        #  [0, 6, 0], [5, 6, 2], [0, 7, 2], [7, 8, 0], [7, 9, 0], [8, 9, 0]],
        # [[1, 2, 0], [1, 3, 0], [1, 4, 0], [1, 5, 0], [1, 6, 0], [1, 8, 0],
        #  [2, 3, 0], [2, 4, 0], [2, 5, 0], [2, 6, 0], [2, 8, 0],
        #  [3, 4, 0], [3, 5, 0], [3, 6, 0],
        #  [4, 5, 0], [4, 6, 0]],
        # [[2, 3, 0], [2, 4, 0], [2, 5, 0], [2, 6, 0], [2, 8, 0],
        #  [3, 4, 0], [3, 5, 0], [3, 6, 0],
        #  [4, 5, 0], [4, 6, 0]],
        # [[0, 1, 0], [0, 2, 0], [0, 3, 0], [1, 2, 0], [1, 3, 0], [2, 3, 0]],
        # [[2, 3, 0], [2, 4, 0], [2, 6, 0],
        #  [3, 4, 0], [3, 6, 0],
        #  [4, 6, 0]]

        # [[0, 1, i], [1, 2, i], [2, 0, j]],
        # [[0, 1, i], [0, 2, j], [1, 2, j]],
        # [[0, 1, j], [0, 2, j], [1, 2, j]],
        # [[0, 1, 1], [0, 2, 0], [0, 3, 0], [1, 2, 0], [1, 3, 0], [2, 3, 2]],
        # [['a', 'c', 1], ['a', 'b', 0], ['a', 'd', 0], ['c', 'b', 0], ['c', 'd', 0], ['b', 'd', 2]],
        # [[0, 1, 0], [0, 2, 0], [0, 3, 0], [0, 4, 0], [1, 4, 2], [2, 3, 0], [2, 4, 2],
        #  [3, 4, 2]],
        # [[0, 1, 0], [0, 2, 0], [0, 3, 0], [0, 4, 0], [0, 5, 0], [0, 6, 0], [0, 7, 0],
        #  [1, 2, 0], [1, 3, 0], [1, 4, 0], [1, 5, 0], [1, 6, 0], [1, 7, 0],
        #  [2, 3, 2], [3, 5, 0], [3, 6, 0], [3, 7, 0], [4, 3, 2], [5, 6, 2], [7, 6, 2]],
        # [[1, 2, 0], [1, 3, 0], [1, 4, 0], [1, 5, 0], [1, 6, 0], [1, 7, 0], [1, 8, 0],
        #  [1, 9, 0],
        #   [2, 3, 0], [2, 4, 0], [2, 5, 0], [2, 6, 0], [2, 8, 0], [2, 7, 2],
        #   [3, 4, 0], [3, 5, 0], [3, 6, 0], [3, 7, 2],
        #   [4, 5, 0], [4, 6, 0],
        #   [7, 8, 0], [7, 9, 2]]
        [[0, 1, 2], [1, 2, 2]], [[0, 1, 2], [0, 2, 2]],
         ]

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

    print('\nGspan canonical graph')
    for g in graphset:
        # g = graphset[4]
        print('\n\tinput graph', g)
        gspan = Gspan(graph=g)
        for _ in range(30):
            map = gspan.graph_randomize()
            gspan.graph_normalize()
            # print('    renormalized graph: {}'.format(gspan.graph))
            glen = len(gspan.graph)
            gspan.minDFS()
            print('\trandomized{}\tminDFS {}'.format(gspan.graph, gspan.mindfs))

    exit(0)
