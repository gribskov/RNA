import sys
import copy


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

    def __init__(self):
        """-----------------------------------------------------------------------------------------
        gspan constructor
        -----------------------------------------------------------------------------------------"""
        self.graph = None
        self.map = None
        self.mindfs = []
        self.g2d = None
        self.d2g = None
        self.unexplored = []

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
        self.d2g = [None for i in range(0, self.vnum)]
        self.g2d = [None for i in range(0, self.vnum)]

        return self.vnum


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

    gspan = Gspan()

    # graph normalization crate an unnormalized graph by doubling the vertex numbers
    print('Graph normalization')
    g = copy.deepcopy(graphset[1])
    print('    original graph: {}'.format(g))
    for edge in g:
        for i in range(0, 2):
            edge[i] *= 2
    print('    un-normalized graph: {}'.format(g))

    gspan.graph = g
    gspan.graph_normalize()
    print('    renormalized graph: {}'.format(gspan.graph))
    if g == gspan.graph:
        print('    passes test')

    exit(0)
