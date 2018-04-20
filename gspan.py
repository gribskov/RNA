import sys
import copy


class Edge(list):
    """=============================================================================================
    Edge class is necessary to implement lexicographic sorting
    e1 < e2 if i1==i2 and j1<j2
               i1<j1 and j1=i2
               e1 <(E,T) and e2 <(E,T) e3 then e1 < e3
    ============================================================================================="""

    def __init__(self, edge=[None, None, None]):
        """-----------------------------------------------------------------------------------------
        An edge is list of vertex 0, vertex 1, and edge type; all ints.  No assumption is made that
        v0 is less than v1.
        -----------------------------------------------------------------------------------------"""
        super(Edge, self).__init__(edge)
        # self = [None, None, None]
        # x=0

    def __eq__(self, other):
        return self[0] == other[0] and self[1] == other[1] and self[2] == other[2]

    def set(self, v0=None, v1=None, e=None):
        """-----------------------------------------------------------------------------------------
        set v0, v1, edge
        :return: True
        -----------------------------------------------------------------------------------------"""
        self[0] = v0
        self[1] = v1
        self[2] = e

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
        self.graph = graph
        self.map = None
        self.vnum = 0
        self.mindfs = []
        self.g2d = None
        self.d2g = None
        self.unexplored = []

        if graph:
            self.graph_normalize()

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

        # convert edge numbers
        for edge in self.graph:
            for i in range(0, 2):
                edge[i] = v.index(edge[i])

        self.vnum = len(v)
        self.map = v

        # initialize d2g and g2d
        self.d2g = [None for _ in range(0, self.vnum)]
        self.g2d = [None for _ in range(0, self.vnum)]

        return self.vnum

    def sort_by_edge(self):
        """-----------------------------------------------------------------------------------------
        sort graph by edges, using current g2d map
        :return:
        -----------------------------------------------------------------------------------------"""
        g2d = self.g2d
        off0 = 256
        off1 = off0 * 256

        def lex(e):
            """
            lexicographic order for edges:
                (known edges not on rightmost path)
                backward from rightmost vertex,
                forward from rightmost vertex,
                forward from internal vertex on rightmost path
            """
            x = 1
            if g2d[e[0]] is not None and g2d[e[1]] is not None:
                # both vertices defined, this is either a known edge or a backward edge
                val = g2d[e[0]] + g2d[e[1]] * off0 + e[2] * off1
                print(e, 'b val', val)
                return val

            elif g2d[e[0]] is not None:
                # e[0] is known, forward extension on rightmost path
                val = g2d[e[0]] + e[2] * off1
                print(e, '0 val', val)
                return val

            elif g2d[e[1]] is not None:
                # e[1] is known, forward extension on rightmost path
                val = g2d[e[1]] * off0 + e[2] * off1
                print(e, '1 val', val)
                return val

            # fallthrough: g2d[e[0]] and g2d[e[1]] are both None:
            # sort by edge type only
            print(e, 'n val', e[2])
            return e[2]

        self.graph.sort(key=lex, reverse=False)

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
                [[0, 1, j], [0, 2, j], [1, 2, j]]]

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
    e.set(1, 2)

    for g in graphset:
        gspan = Gspan(graph=g)
        print(g)
        gspan.g2d = [0, 1, 2]
        gspan.sort_by_edge()
        print(gspan.graph2dfs(), '\n')

exit(0)
