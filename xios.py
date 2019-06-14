class XiosEdge(list):
    """=============================================================================================
    Edge class

    ============================================================================================="""
    # class variables for master definition of edge types
    # a2i   dict    convert edge type alpha to int encoding
    # i2a   list    convert edge type int to alpha encoding
    # arev  dict    reverse edge type in alpha encoding
    # irev  list    reverse edge type in int encoding
    a2i = {'i': 0, 'j': 1, 'o': 2, 's': 3, 'x': 4}
    arev = {'i': 'j', 'j': 'i', 'o': 'o', 's': 3, 'x': 'x'}
    i2a = ['i', 'j', 'o', 's', 'x']
    irev = [1, 0, 2, 3, 4]

    def __init__(self, edge=None):
        """-----------------------------------------------------------------------------------------
        an edge is a triple of  [v0, v1, e]
        v0 originating vertex
        v1 terminating vertex

        -----------------------------------------------------------------------------------------"""
        if edge is None:
            super(XiosEdge, self).__init__([None, None, None])
        else:
            super(XiosEdge, self).__init__(edge)

    def reverse(self):
        """-----------------------------------------------------------------------------------------
        reverse the direction of the edge. Completely agnostic about whether the edges are in alpha
        or integer encoding.  If you are sure the edge type is integer encoded, use flip().

        :return: True
        -----------------------------------------------------------------------------------------"""
        self[0], self[1] = self[1], self[0]
        if isinstance(self[2], int):
            self[2] = XiosEdge.irev[self[2]]
        else:
            self[2] = XiosEdge.arev[self[2]]

        return True

    def flip(self):
        """-----------------------------------------------------------------------------------------
        reverse the direction of the edge. Assumes that the edge type is integer encoded.  Saves
        checking whether the edge type is alpha or int.  If the edge type is not clearly integer,
        use reverse()

        :return: True
        -----------------------------------------------------------------------------------------"""
        self[0], self[1], self[2] = self[1], self[0], XiosEdge.irev[self[2]]

        return True

class Xios(list):
    """=============================================================================================
    A XIOS graph is a description of an RNA topology.
    The basic structure is a list of edges where each edge is [v0, v1, e]
        v0 originating vertex
        v1 terminating vertex
        e edge type (x, i, (j), o, s)
            i edges are directed. i and j indicate forward and backward directions
    
    ============================================================================================="""

    #  __init__ is inherited from list

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


"""
# Encode and look up a DFS code represented as an array or DFS rows, where each row is an
# array of (v1,v2,e).  Each row is encoded  in eight bits, and the entire DFS
# encoded as a hexadecimanl string.  Can only be used with graphs that have seven
# or fewer vertices.
#
# Encoding
#   bit 0 - 1   edge type, i=0, j=1, o=2, x=3
#   bit 2 - 4   to vertex number (0-7)
#   bit 5 - 8   from vertex number (0-7)
#
# USAGE
#    $hex_string = $motif->encodeDfsRowArray( @@dfs );
#------------------------------------------------------------------------------
sub encodeDfsRowArray{
    my ( $motif, @@dfs ) = @@_;
    my $hexstring = "";

    foreach my $dfsrow ( @@dfs ) {
        my $byte;
        my ( $v1, $v2, $e ) = @@{$dfsrow};
        $byte += $v1 << 5;
        $byte += $v2 << 2;
        $byte += $e;

        $hexstring .= sprintf "%02x", $byte;
    }

    return $hexstring;
}

# End of encodeDfsRowArray

"""

if __name__ == '__main__':
    # test that edges are behave properly

    a = ['a', 'b', 'i']
    b = XiosEdge(a)
    print(b)
    b.reverse()
    print('reversed', b)
    b.reverse()
    print('dereversed', b)

    i = [1, 0, 0]
    b = XiosEdge(i)
    print(b)
    b.reverse()
    print('reversed', b)
    b.reverse()
    print('dereversed', b)
    b.flip()
    print('flipped', b)
    b.flip()
    print('deflipped', b)

    c = b
    d = b[:]
    b[1] = 2

    print('a=', a, ' b=', b, ' c=', c, ' d=', d)

    # XIOS
    x = Xios()
    x.append(XiosEdge([0, 1, 0]))
    x.append(XiosEdge([1, 2, 0]))
    x.append(XiosEdge([2, 0, 1]))

