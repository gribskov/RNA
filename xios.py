import sys
import re
import hashlib


class XiosEdge(list):
    """=============================================================================================
    XIOS Edge class
    and edge is a 3-tuple of (v0, v1, edge)
    v0 and v1 can be integers or strings, and need not be consecutive.  edge can be one of the
    allowed characters [ijosx] or the corresponding integer. the authoritative converstions are
    given by the class variables below: a2i, arev, a2a, irev
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

    def copy(self):
        """-----------------------------------------------------------------------------------------
        Creates and returns an independent copy of type XiosEdge

        :return:
        -----------------------------------------------------------------------------------------"""
        return XiosEdge(edge=self[0:3])

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


####################################################################################################
####################################################################################################

class Xios(list):
    """=============================================================================================
    A XIOS graph is a description of an RNA topology.
    The basic structure is a list of edges where each edge is a XIOSEdge, [v0, v1, e]

    XIOSEdge does not constrain whether v0, v1, and e are integers or strings, although e must be an
    integer [0-4] or character[ijosx]

    Most of these variation are fixed up by the normalize method.  the encode/decode methods expect
    the edges to be in all integer form
    ============================================================================================="""

    def __init__(self, **kwargs):
        """-----------------------------------------------------------------------------------------
        Xios graph is a list of XiosEdge

        -----------------------------------------------------------------------------------------"""
        super().__init__()

        for arg in kwargs:
            if arg == 'pair':
                self.len = self.from_pair(kwargs['pair'])
            elif arg == 'PairRNA':
                self.len = self.from_pair(kwargs['PairRNA'])
            elif arg == 'list':
                self.len = self.from_list(kwargs['list'])
            elif arg == 'string':
                self.len = self.from_string(kwargs['string'])
            elif arg == 'graph':
                self.len = self.from_graph(kwargs['graph'])
            elif arg == 'serial':
                self.len = self.from_serial(kwargs['serial'])
            else:
                sys.stderr.write('Xios::__init__ - unknown keyword argument {}'.format(arg))

    def from_list(self, graph):
        """-----------------------------------------------------------------------------------------
        load a graph in the form of a list, into a graph object composed of Edge(). Examples
        [ [0,1,0], [1,2,0], [2,0,1] ]
        [ ['a', 1, 'i'], [1, 2, 'j'], [2, 0, 'j'] ]

        :param graph: list of lists:
        :return: int, number of vertices
        -----------------------------------------------------------------------------------------"""
        for edge in graph:
            # print('type', type(edge[2]))
            if isinstance(edge[2], str):
                edge[2] = XiosEdge.a2i[edge[2]]
            edge = XiosEdge(edge)
            self.append(edge)

        return len(self)

    def from_string(self, graphstr):
        """-----------------------------------------------------------------------------------------
        Reads a graph from a string and converts to non_normalized integers.  all non alphanumerica
        characters are converted to spaces and the string is split.  Values are then taken as
        triples to make the XIOS graph. square brackets and commas are ignored so a python list of
        lists is OK.

        :param graphstr:
        :return: int, number of edges
        -----------------------------------------------------------------------------------------"""
        etype = {'i': 0, 'j': 1, 'o': 2, 's': 3, 'x': 4}

        simple = ''.join([c if c.isalnum() else " " for c in graphstr])
        simple = simple.strip()
        glist = simple.split()

        pos = 0
        triple = []
        for i in glist:
            pos += 1
            if i.isdigit():
                i = int(i)
            else:
                i = etype[i]

            triple.append(i)
            if not pos % 3:
                edge = XiosEdge(triple)
                self.append(edge)
                triple = []

        return len(self)

    def from_pair(self, graph, sedge=False):
        """-----------------------------------------------------------------------------------------
        Converts a PairRNA object to Xios. This requires determining the nesting relationships
        between each pair of stems.

        :param graph: Graph object from RNA/graph.py
        :return: int, number of edges
        -----------------------------------------------------------------------------------------"""
        from topology import PairRNA

        self.clear()

        begin = 0
        end = 1
        for si in range(len(graph)):
            s1 = [graph.pairs[si * 2], graph.pairs[si * 2 + 1]]
            for sj in range(si + 1, len(graph)):
                s2 = [graph.pairs[sj * 2], graph.pairs[sj * 2 + 1]]

                # assuming that s1[begin] is always < s2[begin], i.e., that the PairRNA is canonical

                if s1[end] < s2[begin]:
                    # serial edge
                    if sedge:
                        self.append(XiosEdge([si, sj, 3]))

                elif s1[end] > s2[end]:
                    # i edge, s2 is nested inside s1
                    self.append(XiosEdge([si, sj, 0]))

                else:
                    # pseudoknot
                    self.append(XiosEdge([si, sj, 2]))

        return len(self)

    # def from_pair(self, pair):
    #     """-----------------------------------------------------------------------------------------
    #     creates a Xios object from a PairRNA object, which is defined in topology.py.  In PairRNA
    #     format, each element in the pairlist is a stem, and the pair of values indicate the
    #     locations of the half-stems in the overall structure. For example, [[0,2], [1,3]] is a
    #     pseudoknot.
    #
    #     :param graph:
    #     :return: Xios object
    #     -----------------------------------------------------------------------------------------"""
    #     # convert to pair form
    #     self.from_serial(pair.to_SerialRNA())
    #
    #     return self

    def from_serial(self, graph):
        """-----------------------------------------------------------------------------------------
        creates a Xios graph from a SerialRNA.  SerialRNA is defined in topology.py.  In serial
        format, the positions in a list indicate the left to right positions in the RNA and the
        values indicate which string begins or ends at that position, e.g., [0,1,0,1] for a
        pseudoknot.

        nesting of edges must be determined

        :param graph:
        :return: Xios graph
        -----------------------------------------------------------------------------------------"""
        # convert to pair form
        pos = 0
        pair = [[None, None] for _ in range(len(graph) // 2)]

        for s in graph:
            if pair[s][0] is None:
                pair[s][0] = pos
            else:
                pair[s][1] = pos

            pos += 1

        self.from_graph(pair)
        return self

    def normalize(self):
        """-----------------------------------------------------------------------------------------
        renumber the vertices to be sequential integers.  Convert edge types to integer encoding

        :return: list, map[i] is original vertex name
        -----------------------------------------------------------------------------------------"""
        etype = {'i': 0, 'j': 1, 'o': 2, 's': 3, 'x': 4}

        map = []
        for e in range(len(self)):
            edge = self[e]
            for v in (0, 1):
                if edge[v] not in map:
                    map.append(edge[v])

                edge[v] = map.index(edge[v])

            if type(edge[2]) == 'str':
                if edge[2].isdigit():
                    edge[2] = int(edge[2])
                else:
                    try:
                        edge[2] = etype[edge[2]]
                    except KeyError as err:
                        print(edge, err)

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
            dec2 = int(hex[i + 2:i + 4], 16)
            v0 = (dec1 & 127)
            e = (dec1 & 128) >> 6
            v1 = (dec2 & 127)
            e += (dec2 & 128) >> 7

            # print('[{}, {}, {}]'.format(v0, v1, e))
            edge = XiosEdge([v0, v1, e])
            self.append(edge)

            return True

    def human_encode(self, width=1):
        """-----------------------------------------------------------------------------------------
        A human readable compact encoding for a dfs row. Adaptable to various limits on by changing
        the field width of the vertex integers.

        a row such as 0 1 0 or 0 1 i becomes 0i1 or 00i01 or 000i001 etc.

        :param width: field width for vertices
        :return: str
        -----------------------------------------------------------------------------------------"""
        fmt = '{{:0{}d}}'.format(width)

        string = ''
        for edge in self:

            v0 = fmt.format(edge[0])
            v1 = fmt.format(edge[1])
            e = edge[2]
            if isinstance(edge[2], int):
                e = XiosEdge.i2a[edge[2]]
            string += '{}{}{}.'.format(v0, e, v1)

        return string

    def human_decode(self, code):
        """-----------------------------------------------------------------------------------------
        Decode the encoded graph produced by human_encode(). The current graph is overwritten
        by the new graph coming from the code so this is like reading in a human encoded graph.

        :return: int, number of rows
        -----------------------------------------------------------------------------------------"""
        rowre = re.compile(r'(\d+)([ijosx])(\d+)')
        etype = {'i': 0, 'j': 1, 'o': 2, 's': 3, 'x': 4}

        self.clear()
        for row in code.rstrip('.').split('.'):
            m = rowre.match(row)
            if m:
                (v0, e, v1) = m.group(1, 2, 3)
                e = etype[e]

                edge = XiosEdge([int(v0), int(v1), e])
                self.append(edge)

    def ascii_encode(self):
        """-----------------------------------------------------------------------------------------
        encode a dfs row as two ascii characters, similar to sequence quality in fastq encoding.
        suitable characters are in range 33 (!=21) to 126 (~=7E).  Reserving 4 bits for edge types
        :return:
        -----------------------------------------------------------------------------------------"""
        string = ''
        for edge in self:
            byte1 = 0
            byte1 += (edge[2] & 2) << 5  # edge high order bit (o or s edge)
            byte1 += ((edge[0] + 33) & 63)

            byte2 = 0
            byte2 += (edge[2] & 1) << 6  # low order bit (j or s edge)
            byte2 += (edge[1] + 33 & 63)

            string += '{}{}'.format(chr(byte1), chr(byte2))

        return string

    def ascii_decode(self, code):
        """-----------------------------------------------------------------------------------------
        Decode the ascii encoded graph produced by ascii_encode(). The current graph is
        overwritten by the new graph coming from the code so this is like reading in a ascii
        encoded graph.

        :return: int, number of rows
        -----------------------------------------------------------------------------------------"""
        self.clear()
        n = len(code)
        for i in range(0, len(code), 2):
            dec1 = ord(code[i])
            dec2 = ord(code[i + 1])
            v0 = (dec1 & 63) - 33
            e = (dec1 & 64) >> 5
            v1 = (dec2 & 63) - 33
            e += (dec2 & 64) >> 6

            # print('[{}, {}, {}]'.format(v0, v1, e))
            edge = XiosEdge([v0, v1, e])
            self.append(edge)

        return True


import sys
import copy
import random
from functools import total_ordering
from xios import XiosEdge, Xios

"""=================================================================================================
MotifDB class is for creating and using XIOS graph dictionaries


Michael Gribskov     15 June 2019
================================================================================================="""
import json
import datetime
from xios import Xios


class MotifDB():
    """=============================================================================================
    I keep going back and forth over whether the object should be a dict, or whether it should have
    defined fields.  Defined fields requires fewer accessors, but a general accessor with getattr
    could be used.
        motifdb.fields = list of fields in order for serialization
        motifdb.information = metadata, standard info has accessors, but any notes can be added
        motifdb.db = {dfs_string:n_stem, ...}
        MotifDB.parent = {dfsstring:[parent_dfsstring0-, parent_dfsstring1, ...], ...

    The motif database is indexed by the minimum dfs string. The dfs string could be text (current)
    or use previous hexadecimal encodings to save space. Since the minimum dfs is unique, a
    dictionary can be used.
    ============================================================================================="""

    def __init__(self, **kwds):
        self.fields = ['information', 'db', 'parent']
        self.information = {}  # for metadata
        self.db = {}
        self.parent = {}

        for key in kwds:
            if key == 'json':
                self.fromJSON(kwds[key])
            else:
                sys.stderr.write('MotifDB::init - unknown keyword ({})'.format(key))

    def add_with_len(self, motif, n_stems):
        """-----------------------------------------------------------------------------------------
        Add a single motif to the database

        :param motif: string, text string describing the dfs
        :param n_stems: int, number of stems (not edges) in the motif
        :return: int, number of motifs in database
        -----------------------------------------------------------------------------------------"""
        self.db[motif] = n_stems
        self.n = len(self.db)

        return self.n

    def pickle(self, fh):
        """-----------------------------------------------------------------------------------------
        Pickle the object and save to file.  Use filehandle instead of filename so that multiple
        objects can be pickled

        :param fh: filehandle, open file for writing
        :return: True
        -----------------------------------------------------------------------------------------"""
        import pickle

        pickle.dump(self, fh)

        return

    @staticmethod
    def unpickle(filename):
        """-----------------------------------------------------------------------------------------
        Restore pickle database

        :param filename:
        :return: boolean, True if sucessful
        -----------------------------------------------------------------------------------------"""
        import pickle

        fh = None
        try:
            fh = open(filename, 'rb')
        except OSError:
            sys.stderr.write(f"MotifDB - can't open pickle file ({filename})")

        return pickle.load(fh)

    def setdate(self):
        """-----------------------------------------------------------------------------------------
        set date in information
        :return: datetime
        -----------------------------------------------------------------------------------------"""
        daytime = datetime.datetime.now()
        self.information['date'] = daytime.strftime('%Y-%m-%d %H:%M:%S')

        return daytime

    def setname(self, string):
        """-----------------------------------------------------------------------------------------
        set database name in information

        :return: int, length of new name field
        -----------------------------------------------------------------------------------------"""
        self.information['name'] = string

        return len(string)

    def setsource(self, source):
        """-----------------------------------------------------------------------------------------
        Name of the program creating the database
        :return:
        -----------------------------------------------------------------------------------------"""
        self.information['source'] = source

        return len(source)

    def checksum(self):
        """-----------------------------------------------------------------------------------------
        calculate a checksum.  The checksum is based on the JSON string of the db and parent
        so it should not be affected by changes to comments.

        :return: str, hexadecimal md5 checksum
        -----------------------------------------------------------------------------------------"""
        # content = json.dumps(self.db) + json.dumps(self.parent)
        # result = hashlib.md5(content.encode())
        from functools import reduce
        import zlib
        result = reduce(lambda x, y: x ^ y, [zlib.adler32(bytes(repr(t), 'utf-8')) for t in self.db.items()])
        result = result ^ reduce(lambda x, y: x ^ y, [zlib.adler32(bytes(repr(t), 'utf-8')) for t in
                                                      self.parent.items()])

        return result
    # try
    #     return hashlib.md5(content.encode()).hexdigest()
    # or sha256

    def toJSON(self):
        """-----------------------------------------------------------------------------------------
        Convert database to JSON string
        self.fields specifies the order of the fields in the output

        :return: str
        -----------------------------------------------------------------------------------------"""
        fields = self.fields

        self.information['checksum'] = self.checksum()
        dispatch = []
        for i in range(len(fields)):
            # look up the values of each field and add to dispatch list
            dispatch.append(getattr(self, fields[i]))

        return json.dumps({fields[i]: dispatch[i] for i in range(len(fields))}, indent=4)

    def fromJSON(self, fpin):
        """-----------------------------------------------------------------------------------------
        Read in a JSON serialized MotifDB written by toJSON

        :param fpin:  file/str, either a file open for reading, or a string giving the path
        :return: True
        -----------------------------------------------------------------------------------------"""
        # check the argument, if it is a file pointer (open file) it can just be used
        # otherwise, open the file
        fp = None
        if isinstance(fpin, str):
            try:
                fp = open(fpin, 'r')
            except OSError:
                sys.stderr.write('MotifDB.fromJSON - error opening file ({})'.format(fpin))
        else:
            fp = fpin

        # read the JSON
        j = json.load(fp)

        for i in self.fields:
            # look up the values of each field and add to dispatch list
            setattr(self, i, j[i])

        old_checksum = self.information['checksum']
        self.information['checksum'] = self.checksum()
        checksum_is_ok = old_checksum == self.information['checksum']
        if checksum_is_ok is False:
            sys.stderr.write('MotifDB::fromJSON - checksums do not match (file:{} new:{}'.
                             format(old_checksum, self.information['checksum']))

        return checksum_is_ok

    def toFile(self, fp):
        """-----------------------------------------------------------------------------------------
        Write a formatted version of the motif database to a file. use sys.stdout as the file if
        you want it on STDOUT.  The JSON format is authoritative, this is for reference only.

        :param fp:
        :return:
        -----------------------------------------------------------------------------------------"""
        fields = self.fields

        self.information['checksum'] = self.checksum()

        for i in range(len(fields)):
            data = getattr(self, fields[i])
            if fields[i] == 'information':
                # metadata fields
                for tag in data:
                    fp.write(f'{tag}\t{data[tag]}\n')
                # number of motifs
                fp.write(f'motifs:\t{len(self.db)}\n')
            elif fields[i] == 'db':
                for m in data:
                    fp.write(f'\t{self.db[m]}\t{m}\n')
                    for p in self.parent[m]:
                        fp.write(f'\t\t\t{p}\n')

        return


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
        reverse the direction of the edge, edge type 0 (i) and 1 (j) are inverses of each other.
        o and s edges are undirected so just vertex 0 and 1 change.

        :return: True
        -----------------------------------------------------------------------------------------"""
        # swap vertex 0 and 1
        self[0], self[1] = self[1], self[0]
        if self[2] < 2:
            # ^ is XOR
            self[2] ^= 1

        return True

    def copy(self):
        return Edge(self)


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

        unexplored, the stack of partial solutions stores [g2d, edge, row], see save()/restore()
        -----------------------------------------------------------------------------------------"""
        from topology import PairRNA

        self.graph = Xios()  # graph in normalized labelling
        self.nforward = 0  # number of forward edges in sorted graph
        self.nbackward = 0  # number of backward edges in sorted graph
        self.nunknown = 0  # number of unknown edges in sorted graph
        self.vnum = 0  # number of vertices in graph
        self.vnext = 0  # next dfs vertex available to use
        self.row = 0
        self.mindfs = Xios()
        self.mindfslen = 0
        self.mindfsg2d = []
        self.g2G = []  # list to convert g labels (indices) to original labels
        self.g2d = []  # list to covert g labels to dfs labels
        self.d2g = []  # list to conver dfs labels to g labels
        self.unexplored = []  # stack of partial solutions that need to be searched

        if graph:
            if isinstance(graph, PairRNA):
                self.from_pair(graph)

            elif isinstance(graph, list):
                self.from_list(graph)

            elif isinstance(graph, str):
                self.from_string(graph)

            elif isinstance(graph, Xios):
                self.graph = graph

            else:
                sys.stderr.write('Gspan::__init__ - unknown graph type ({})\n'.format(graph))

            self.graph_normalize()

    def from_list(self, graph):
        """-----------------------------------------------------------------------------------------
        load a graph in the form of a list, into a graph object composed of Edge()
        TODO use Xios object

        :param graph: list of lists:
        :return: int, number of vertices
        -----------------------------------------------------------------------------------------"""
        self.graph = list()
        v = 0
        for edge_in in graph:
            edge = Edge(edge_in)
            self.graph.append(edge)
            v += 1

        for _ in self.graph:
            self.mindfs.append(XiosEdge())

        self.vnum = v

        return v

    def from_string(self, graphstr):
        """-----------------------------------------------------------------------------------------
        Reads a graph from a string and converts to non_normalized integers.  For now graphs are
        a list of edges, each edge is a triple of vertex0 vertex1 edge_type.  All are abritrary
        strings without quotation marks. square brackets and commas are ignored so a python list of
        lists is OK

        TODO use Xios object

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

    def from_pair(self, pair_rna):
        """-----------------------------------------------------------------------------------------
        Read a topology as PairRNA object (see Topology). since it is stored in Gspan as a Xios
        object, read the PairRNA object into Xios and store in Gspan

        :param pair: PairRNA
        :return: True
        -----------------------------------------------------------------------------------------"""
        xios = Xios(pair=pair_rna)
        self.graph = xios

        return True

    def flip(self, row=0):
        """-----------------------------------------------------------------------------------------
        convert all j edges to i edges from row to end of graph. Not currently used

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

        TODO move to Xios object

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
            for i in range(2):
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
        randomly relabel the graph vertices.  this is useful for generating graphs that have the
        same canonical form.

        TODO move to Xios object

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

    def save(self, edge):
        """-----------------------------------------------------------------------------------------
        push a partial solution onto the unexplored stack. Stores copy of edge, row

        :return: integer, length of stack
        -----------------------------------------------------------------------------------------"""
        self.unexplored.append([edge, self.row])

        return len(self.unexplored)

    def restore(self):
        """-----------------------------------------------------------------------------------------
        pop a saved dfs code from the unexplored stack
        swap the edge with the current edge at the specified row

        :return: integer, next available dfs vertex
        -----------------------------------------------------------------------------------------"""
        while len(self.unexplored) > 0:
            edge, self.row = self.unexplored.pop()

            # this makes the popped edge current edge by moving it to the restored row
            try:
                # if the popped edge can't be found it must be in the reverse orientation
                epos = self.graph.index(edge)
            except ValueError:
                t = edge[2]
                if t < 2:
                    t ^= 1
                epos = self.graph.index([edge[1], edge[0], t])

            # move the popped edge to the desired row in the graph
            self.graph[epos] = self.graph[self.row]
            self.graph[self.row] = edge
            self.row += 1

            # rebuild g2d and d2g from scratch
            self.vnext = 0
            self.d2g = [None for _ in self.d2g]
            self.g2d = [None for _ in self.g2d]
            for i in range(self.row):
                edge = self.graph[i]
                for e in range(0, 2):
                    if edge[e] not in self.d2g:
                        self.g2d[edge[e]] = self.vnext
                        self.d2g[self.vnext] = edge[e]
                        self.vnext += 1

            # check to see if the graph could be minimum, if not, go on to the next one
            if not self.minimum(0):
                continue

            # restored graph could be minimum
            return True

        # end of loop over unexplored stack, if you reach here, the stack is empty
        return False

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

    def graph2dfs(self, first):
        """-----------------------------------------------------------------------------------------
        convert the node labels in the graph to dfs labels and store in mindfs

        :first: integer, first row to convert/copy
        :return: list of edges
        -----------------------------------------------------------------------------------------"""
        row = self.row
        g2d = self.g2d
        dfs = self.mindfs
        for i in range(first, row):
            edge = self.graph[i]
            try:
                dfs[i] = [g2d[edge[0]], g2d[edge[1]], edge[2]]
            except IndexError:
                # dfs needs to be extended
                dfs.append([g2d[edge[0]], g2d[edge[1]], edge[2]])

        self.mindfslen = self.row

        return dfs

    def initDFS(self):
        """-----------------------------------------------------------------------------------------
        stores all possible initial edges to unexplored list

        :return: int, size of unexplored list
        -----------------------------------------------------------------------------------------"""
        self.row = row = 0
        self.sort(begin=row)
        self.unexplored = []

        first_edge_type = self.graph[row][2]
        if first_edge_type == 2:
            # special case: when all edges are undirected, save both orientations on unexplored
            for edge in self.graph:
                if edge[2] != first_edge_type:
                    break
                self.save(edge.copy())
                # erev = self.graph[row].copy()
                erev = self.graph[row].copy()
                erev.reverse()

                self.save(erev)
        else:
            # first edge type is not undirected, save just one orientation
            for edge in self.graph:
                if edge[2] != first_edge_type:
                    break
                self.save(edge.copy())

        return len(self.unexplored)

    def minDFS(self):
        """-----------------------------------------------------------------------------------------
        Find the minimum DFS code by exhaustive search from all possible initial edges.

        :return:
        -----------------------------------------------------------------------------------------"""
        if len(self.graph) < 1:
            return []

        self.initDFS()
        searching = self.restore()
        graph = self.graph

        while searching:
            # sort
            first = self.row
            self.sort(begin=self.row)

            # add backward edges, since the graph is sorted, all this requires is updating the
            # begin point in the graph array. backward edges never add a new vertex so d2g and g2d
            # are untouched
            self.row += self.nbackward

            if self.row < len(graph):
                # skip adding forward edges if done

                # save equivalent forward edges.  We know the number of forward edges from the sort,
                # the edges are equivalent if they have the same v0 and edge type (v2)
                first_edge_type = graph[self.row][2]
                v0 = self.g2d[graph[self.row][0]]
                for edge in graph[self.row + 1:]:
                    if edge[2] != first_edge_type:
                        break
                    if self.g2d[edge[0]] != v0:
                        break

                    # if you pass these tests, fall though to saving this edge on stack
                    self.save(edge.copy())

                # add the forward extension to g2d and d2g
                self.d2g[self.vnext] = graph[self.row][1]
                self.g2d[graph[self.row][1]] = self.vnext
                self.vnext += 1
                self.row += 1

            # check to see if the graph is a possible minimum, or if this graph is done
            if not self.minimum(first) or self.row == len(graph):
                searching = self.restore()

        return self.mindfs

    def minimum(self, first):
        """-----------------------------------------------------------------------------------------
        check if the current dfs is <= the min dfs (returns True). otherwise, returns False.
        mindfs is stored in d space, current dfs must be compared in d space

        :return: logical
        -----------------------------------------------------------------------------------------"""

        row = self.row
        if row > self.mindfslen:
            self.graph2dfs(first)
            return True

        cmp = None
        for i in range(first, row):
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
                        sys.stdout.write('\t{}\n'.format(self.graph))
                        sys.stdout.write('\t{}\n'.format(self.mindfs))
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
                elif cedge[1] > medge[1]:
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
            self.graph2dfs(first)
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
        try:
            if edge[0] < edge[1]:
                direction = 'f'
        except TypeError:
            print(f'edge0:{edge[0]}   edge1: {edge[1]}')

        return direction

    def edge_g2d(self, edge):
        """-----------------------------------------------------------------------------------------
        convert edge from g space to d space

        :param edge: Edge (g space)
        :return: Edge (d space)
        -----------------------------------------------------------------------------------------"""
        return Edge([self.g2d[edge[0]], self.g2d[edge[1]], edge[2]])


# ==================================================================================================
# testing
# ==================================================================================================


if __name__ == '__main__':

    # the main testing program is at the very end

    def test_XiosEdge():
        """-----------------------------------------------------------------------------------------
        Test SiosEdge class
        :return:
        -----------------------------------------------------------------------------------------"""
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

        e3 = e2.copy()
        e3[0], e3[1] = e3[1], e3[0]
        e3[0] = 5
        print('copy {} | {}'.format(e2, e3))

        if e1 < e2:
            print('e1 smaller')
        if e2 < e1:
            print('e2 smaller')

        return


    def test_Xios():
        """-----------------------------------------------------------------------------------------
        Test Xios class

        :return:
        -----------------------------------------------------------------------------------------"""
        print('Testing Xios class')
        # mapping of edge types to integers to avoid quoting all the time
        i = 0
        j = 1
        o = 2
        s = 3

        print('\nLoad from python list')
        rnas = [[[0, 1, 0], [1, 2, 0], [2, 0, 1]], [['a', 1, 'i'], [1, 2, 'j'], [2, 0, 'j']]]
        for rna in rnas:
            x = Xios(list=rna)
            print('\t{}'.format(x))

        print('\nLoad from python list-like string')
        rnas = ['[[0, 1, 0], [1, 2, 0], [2, 0, 1]]', '(0,1,0) (1,2,0) (2,0,1)', '0 1 i 1 2 i 2 0 j']
        for rna in rnas:
            x = Xios(string=rna)
            print('\t{}\t =>\t{}'.format(rna, x))

        print('\nRead from Pair format (no normalization)')
        rnapairs = [[[0, 1], [2, 3]],
                    [[0, 2], [1, 3]],
                    [[0, 3], [1, 2]],
                    [[0, 5], [1, 2], [3, 4]],
                    [[1, 7], [1, 3], [2, 5], [4, 7]]]
        for rna in rnapairs:
            xg = Xios(graph=rna)
            print('\tinput:{}\tXios:{}'.format(rna, xg))

        print('\nRead from SerialRNA format (no normalization)')
        rnas = [[0, 1, 2, 2, 1, 0], [1, 0, 2, 2, 0, 1], [0, 1, 2, 1, 2, 0], [2, 0, 1, 2, 0, 1]]
        for rna in rnas:
            xg = Xios(serial=rna)
            print('\tinput:{}\tXios:{}'.format(rna, xg))

        print('\nNormalization')

        graph = [['a', 13, 'i'], [11, 12, 'i'], [12, 'a', 'j']]
        print('\tx2 read from list', graph)
        x2 = Xios(list=graph)
        map = x2.normalize()
        print('\tnormalized={}\tmap:{}\n'.format(x2, map))

        graph = '[[0, 1, 0], [1, 2, 0], [2, 0, 1]]'
        print('\tx3 read from string', graph)
        x3 = Xios(string=graph)
        map = x3.normalize()
        print('\tnormalized={}\tmap:{}\n'.format(x3, map))

        # encoding/decoding

        print('\nTest encoding')
        test = Xios()
        x4 = Xios(graph=rnapairs[4])

        graphs = [x2, x3, x4]
        print('\thuman encoding')
        for rna in graphs:
            print('\t\tinput {}'.format(rna))
            code = rna.human_encode()
            print('\t\tcoded {}'.format(code))
            test.human_decode(code)
            print('\t\tdecoded {}\n'.format(test))

        print('\thex2 encoding')
        for rna in graphs:
            print('\t\tinput {}'.format(rna))
            code = rna.hex2_encode()
            print('\t\tcoded {}'.format(code))
            test.hex2_decode(code)
            print('\t\tdecoded {}\n'.format(test))

        print('\tascii encoding')
        for rna in graphs:
            print('\t\tinput {}'.format(rna))
            code = rna.ascii_encode()
            print('\t\tcoded {}'.format(code))
            test.ascii_decode(code)
            print('\t\tdecoded {}\n'.format(test))

        exit(1)

        # graph normalization create an unnormalized graph by doubling the vertex numbers
        # XIOS
        x = Xios()
        x.append(XiosEdge([0, 1, 0]))
        x.append(XiosEdge([1, 2, 0]))
        x.append(XiosEdge([2, 0, 1]))

        return


    def test_Gspan():
        """-----------------------------------------------------------------------------------------
        Test Gspan class
        :return:
        -----------------------------------------------------------------------------------------"""
        '''
         graphset(0:2) have the same canonical representation as graphset[0]
         '''
        graphset = [
            # [[0,1,0], [0,2,0], [1,2,0], [0,3,2], [1,3,2], [2,3,0], [3,4,0]],
            # [[0,1,0], [0,2,2], [0,4,2], [1,2,2], [2,3,2], [2,4,0], [2,5,2], [4,5,2], [5,3,0]],
            # [[0,1,0], [0,2,0], [0,3,0], [0,4,0], [0,5,0], [1,2,0], [4,5,0]],
            # [[0, 1, 0], [0, 2, 0], [0, 3, 0], [0, 4, 0], [0, 5, 0], [1, 2, 0], [4, 5, 0], [0,6,0],
            # [1,6,2], [4,6,2]],
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
            # [[0, 1, 2], [1, 2, 2]], [[0, 1, 2], [0, 2, 2]],
            # [[0, 1, 2], [0, 2, 0], [0, 3, 2], [0, 4, 2], [1, 2, 0], [1, 3, 0], [1, 4, 0], [1, 5, 0],
            #  [3, 4, 2], [4, 5, 2]],
            # [[0, 1, 0], [0, 2, 0], [0, 3, 0], [0, 4, 2], [0, 5, 0], [1, 2, 2], [2, 3, 2], [2, 4, 2],
            #  [3, 4, 2], [4, 5, 0]],
            # [[0, 1, 2], [0, 2, 2], [1, 2, 2], [1, 3, 2], [2, 3, 2], [3, 4, 2], [4, 5, 2]],
            # [[0, 1, 2], [1, 2, 2], [2, 3, 2], [2, 4, 2], [3, 4, 2], [3, 5, 2], [4, 5, 2]],
            # [[0, 1, 2], [1, 2, 2], [2, 3, 2], [2, 4, 2], [3, 4, 2], [4, 5, 2]],
            # [[0, 1, 2], [1, 2, 2], [1, 3, 2], [2, 3, 2], [3, 4, 2], [4, 5, 2]],
            # [[0, 1, 2], [1, 2, 2], [1, 3, 2], [2, 3, 2], [3, 4, 2], [4, 5, 2], [5, 6, 2]],
            # [[0, 1, 2], [1, 2, 2], [2, 3, 2], [3, 4, 2], [3, 5, 2], [4, 5, 2], [5, 6, 2]],
            [[1, 5, 0], [1, 3, 0], [1, 4, 0], [1, 2, 0], [1, 0, 2], [2, 0, 2]]
            ]
        # testing reading graphs from string
        gstr = '[[a, c, 1], [a, b, 0], [a, d, 0], [c, b, 0], [c, d, 0], [b, d, 2]]'
        # g = Gspan(gstr)

        for graph in graphset:
            print('\nGraph normalization')
        g = copy.deepcopy(graph)
        print('    original graph: {}'.format(g))

        # for edge in g:
        #     for i in range(0, 2):
        #         edge[i] *= 2
        # print('    un-normalized graph: {}'.format(g))

        print('\nGspan canonical graph')
        for g in graphset:
            # g = graphset[4]
            print('\n\tinput graph', g)
            gspan = Gspan(graph=g)
            for _ in range(1):
                # map = gspan.graph_randomize()
                gspan.graph_normalize()
                # print('    renormalized graph: {}'.format(gspan.graph))
                glen = len(gspan.graph)
                gspan.minDFS()
                print('\trandomized{}\n\tminDFS {}'.format(gspan.graph, gspan.mindfs))

        return


    def test_MotifDB():
        """-----------------------------------------------------------------------------------------
        testing - motifDB

        :return:
       ----------------------------------------------------------------------------------------- """

        motif = MotifDB()
        motif.setdate()
        motif.setname('test database')

        motif.add_with_len('aabbcc', 1)
        motif.add_with_len('bbccdd', 2)
        motif.add_with_len('ccddee', 2)
        print(motif.toJSON())

        return


    # ##############################################################################################
    # Testing
    # ##############################################################################################

    def test_Main():
        # test_XiosEdge()
        # test_Xios()
        # test_Gspan()
        test_MotifDB()

        return


    test_Main()
    exit(0)
