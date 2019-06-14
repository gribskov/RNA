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

    def int(self):
        """-----------------------------------------------------------------------------------------
        convert the edge type to integer encoding

        :return: True
        -----------------------------------------------------------------------------------------"""
        if not isinstance(self[2], int):
            self[2] = XiosEdge.a2i[self[2]]

        return


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
        load a graph in the form of a list, into a graph object composed of Edge(). Examples
        [ [0,1,0], [1,2,0], [2,0,1] ]
        [ ['a', 1, 'i'], [1, 2, 'j'], [2, 0, 'j'] ]

        :param graph: list of lists:
        :return: int, number of vertices
        -----------------------------------------------------------------------------------------"""
        for edge in graph:
            edge = XiosEdge(edge)
            self.append(edge)
        return len(self)

    def from_string(self, graphstr):
        """-----------------------------------------------------------------------------------------
        Reads a graph from a string and converts to non_normalized integers.  all non alphanumerica
        characters are converted to spaces and the strin is split.  Values are then taken as triples
        to make the XIOS graph. square brackets and commas are ignored so a python list of
        lists is OK

        :param graphstr:
        :return:
        -----------------------------------------------------------------------------------------"""
        simple = ''.join([c if c.isalnum() else " " for c in graphstr])
        simple = simple.strip()
        graph = simple.split()
        for i in range(0, len(graph), 3):
            edge = XiosEdge([graph[1], graph[i + 2], graph[i + 2]])
            self.append(edge)

        return len(self)

    def normalize(self):
        """-----------------------------------------------------------------------------------------
        renumber the vertices to be sequential integers.  Convert edge types to integer encoding

        :return: list, map[i] is original vertex name
        -----------------------------------------------------------------------------------------"""
        n = 0
        map = []
        for edge in self:
            for v in (0, 1):
                if edge[v] not in map:
                    map.append(edge[v])
                    n += 1

                edge[v] = map.index(edge[v])

            if edge[2].isdigit():
                edge[2] = int(edge[2])
            else:
                edge.int()

        return map

    def hex_encode(self):
        """-----------------------------------------------------------------------------------------
        endcode the entire graph as a hexadecimal string.  Each row is encoded  in eight bits, using
        3 bits for each vertex, so it can only be used with graphs that have seven or fewer vertices.

        Encoding
          bit 0 - 1   edge type, i=0, j=1, o=2, x=3
          bit 2 - 4   to vertex number (0-7)
          bit 5 - 8   from vertex number (0-7)

        :return:
        -----------------------------------------------------------------------------------------"""
        hexstring = ''
        for edge in self:
            byte = 0
            byte += edge[0] << 5
            byte += edge[1] << 2
            byte += edge[2]

            hexstring += '{:02x}'.format(byte)

        return hexstring

    def hex_decode(self, hex):
        """-----------------------------------------------------------------------------------------
        Decode the hexadecimal encoded graph produced by hex_encode(). The current graph is
        overwritten by the new graph coming from the hex code so this is like reading in a hex
        encoded graph.

        :return: int, number of rows
        -----------------------------------------------------------------------------------------"""
        self.clear()
        n = 0
        for i in range(0, len(hex), 2):
            n += 1
            dec = int(hex[i:i + 2], 16)
            v0 = (dec & 224) >> 5
            v1 = (dec & 28) >> 2
            e = (dec & 3)

            # print('[{}, {}, {}]'.format(v0, v1, e))
            edge = XiosEdge([v0, v1, e])
            self.append(edge)

        return n

    def hex2_encode(self):
        """-----------------------------------------------------------------------------------------
        endcode the entire graph as a hexadecimal string using two bytes.  Each vertex is encoded as
        7 bits so it can only be used with graphs that have 127 or fewer vertices.

        Encoding
          bit 0       edge type part 1, i=0, j=1, o=2, x=3
          bit 1 - 7   v0 (0-128)
          bit 0       edge type part 2, i=0, j=1, o=2, x=3
          bit 1 - 7   v1 (0-128)


        :return:
        -----------------------------------------------------------------------------------------"""
        hexstring = ''
        for edge in self:
            byte1 = 0
            byte1 += (edge[2] & 2) << 7  # edge high order bit (o or s edge)
            byte1 += (edge[0] & 127)

            byte2 = 0
            byte2 += (edge[2] & 1) << 7  # low order bit (j or s edge)
            byte2 += (edge[1] & 127)

            hexstring += '{:02x}{:02x}'.format(byte1, byte2)

        return hexstring

    def hex2_decode(self, hex):
        """-----------------------------------------------------------------------------------------
        Decode the hexadecimal encoded graph produced by hex2_encode(). The current graph is
        overwritten by the new graph coming from the hex code so this is like reading in a hex
        encoded graph.

        :return: int, number of rows
        -----------------------------------------------------------------------------------------"""
        self.clear()
        n = 0
        for i in range(0, len(hex), 4):
            n += 1
            dec1 = int(hex[i:i + 2], 16)
            dec2 = int(hex[i+2:i + 4], 16)
            v0 = (dec1 & 127)
            e = (dec1 & 128) >> 6
            v1 = (dec2 & 127)
            e += (dec2 & 128) >> 7

            print('[{}, {}, {}]'.format(v0, v1, e))
            edge = XiosEdge([v0, v1, e])
            self.append(edge)



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
    print('x=', x)

    graph = [[0, 1, 0], [1, 2, 0], [2, 0, 1]]
    print('\nx1 read from list', graph)
    x1 = Xios()
    x1.from_list(graph)
    print('x1=', x1)

    graph = [['a', 13, 'i'], [11, 12, 'i'], [12, 'a', 'j']]
    print('\nx2 read from list', graph)
    x2 = Xios()
    x2.from_list(graph)
    print('x2=', x2)

    graph = '[[0, 1, 0], [1, 2, 0], [2, 0, 1]]'
    print('\nx3 read from string', graph)
    x3 = Xios()
    x3.from_string(graph)
    print('x3=', x3)

    print('\nrnormalize x3')
    print('x=', x3)
    map = x3.normalize()
    print('normalized=', x3, map)

    print('\nrnormalize x2')
    print('x=', x2)
    map = x2.normalize()
    print('normalized=', x2, map)

    print('\nhex encoding x1')
    hex = x1.hex2_encode()
    print(x1, hex)
    x1.hex2_decode(hex)
    print('decoded', x1)
    print('\nhex encoding x2')
    hex = x2.hex2_encode()
    print(x2, hex)
    x2.hex2_decode(hex)
    print('decoded', x2)
