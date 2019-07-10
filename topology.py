"""=================================================================================================
    There are several ways to represent a folded RNA structure as a linear string at an abstract
    level where we consider simply paired regions rather than individual base pairs.  These include
    the Giegerich abstract shapes approach, where the runs of parentheses and dots in a Vienna
    formatted structure are condensed to single brackets (square brackets in their formulation).
    This package handles the manipulation of these kinds of linear representations (including some
    extensions for pseudoknots). A second representation numbers the stem regions based on the
    pairing.  Again, this can be condensed so each stem is represented by a single digit
    For instance, for a tRNA

    >S.cerevisiae_tRNA-PHE M10740/1-73
    GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA
    (((((((..((((........)))).((((.........)))).....(((((.......)))))))))))).   Vienna format
    1111111  2222        2222 3333         3333     44444       444441111111    stem number
    1111111  2222        3333 4444         5555      6666       777788888888    position

    abstract shape (level 3) [ [ ] [ ] [ ] ]    (spaces added for readability)
    serial format   1 2 2 3 3 4 4 1             index indicates position, value indicates stem
    pair format   (1,8) (2,3) (4,5) (6,7)       each item is a stem, values are the positions
                  1 8 2 3 4 5 6 7               (serialized pair format)

    example with a pseudoknot added
    GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA
    (((((((..((((........)))).((((.........)))).[[[.(((((..]]]..)))))))))))).   Vienna format
    1111111  2222        2222 3333         3333 444 55555  444  555551111111    stem number
    1111111  2222        3333 4444         5555 666 77777  888  999991111111    position
                                                                     0000000
    ( ( ) ( ) [ ( ] ) )
    1 2 2 3 3 4 5 4 5 1
    (1, 10) (2, 3) (4, 5) (6, 8) (7, 9)
    1 10 2 3 4 5 6 8 7 9

    In addition to these simple formats, there are more complicated outputs written by structure
    prediction programs such as :

    CT file
        mfold and RNA Structure

    RNAml
================================================================================================="""
import sys
import copy
from datetime import datetime
import random
from lxml import etree
from xios import Xios


####################################################################################################
####################################################################################################
class Topology:
    """=============================================================================================
    an RNA topology describes an RNA with one more of the following attributes:
        list of stems (with start, stop, and vienna)
        vienna (dot-bracket) string
        nucleic acid sequence
        sequence ID
        sequence description (documentation)
        sequence length
        comments

    many topologies will have only some of these.

    Many functions of topology.pm in the Perl version have not yet been ported
    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        self.stem_list = []
        self.edge_list = []
        self.adjacency = []
        self.vienna = ''
        self.sequence = ''
        self.sequence_id = ''
        self.sequence_doc = ''
        self.sequence_length = 0
        self.comment = []

    def format_stem_list(self):
        """-----------------------------------------------------------------------------------------
        Construct a formatted string in which the stemlist columns line up nicely. stem_list is a
        list of dict with fields
            name
            center
            left_begin
            left_end
            right_begin
            right_end
            left_vienna
            right_vienna
        
        :return: str
        -----------------------------------------------------------------------------------------"""
        string = ''
        colmax = {}
        for s in self.stem_list:
            for column in s:
                colstr = '{}'.format(s[column])
                if column in colmax:
                    colmax[column] = max(colmax[column], len(colstr))
                else:
                    colmax[column] = len(colstr)

        fmt = ' {{:>{}}}  {{:{}}}  [ {{:{}}} {{:{}}} {{:{}}} {{:{}}} ]  {{:>{}}}  {{:{}}}/n'. \
            format(colmax['name'],
                   colmax['center'],
                   colmax['left_begin'],
                   colmax['left_end'],
                   colmax['right_begin'],
                   colmax['right_end'],
                   colmax['left_vienna'],
                   colmax['right_vienna'])

        for s in self.stem_list:
            string += fmt.format(s['name'],
                                 s['center'],
                                 s['left_begin'],
                                 s['left_end'],
                                 s['right_begin'],
                                 s['right_end'],
                                 s['left_vienna'],
                                 s['right_vienna'])

        return string.rstrip('/n')

    def format_edge_list(self, fieldwidth=4):
        """-----------------------------------------------------------------------------------------
        edge_list is written from adjacency list

        :return: str
        -----------------------------------------------------------------------------------------"""
        string = ''
        adj = self.adjacency
        fmtrow = '{{:>{}}}:'.format(fieldwidth - 1)
        fmtcol = '{{:>{}}}{{}}'.format(fieldwidth - 1)

        r = 0
        for row in adj:
            string += fmtrow.format(r)
            string += ''.join(
                [fmtcol.format(i, row[i]) for i in range(r + 1, len(row)) if row[i] != 's'])
            string += '/n'
            r += 1

        return string.strip('/n')

    def format_adjacency(self, fieldwidth=3):
        """-----------------------------------------------------------------------------------------
        Return a formatted string with the adjacency matrix.

        :param fieldwidth: int, with of field in formatted string
        :return: str
        -----------------------------------------------------------------------------------------"""
        string = ' ' * fieldwidth
        adj = self.adjacency

        fmt = '{{:>{}}}'.format(fieldwidth)
        string += ''.join([fmt.format(i) for i in range(len(adj))])
        r = 0
        for row in adj:
            string += '/n'
            string += fmt.format(r)
            string += ''.join([fmt.format(row[i]) for i in range(len(adj))])
            r += 1

        return string

    @staticmethod
    def iwrite(fp, string, level, tag='', indent_len=2, indent_char=' '):
        """-----------------------------------------------------------------------------------------
        write a string, possibly multi-line, with indentation.  If tag is provided, opening and
        closing XML tags are written before and after the block (on separate lines at the same 
        indentation level).

        :param fp: file, open file for writing
        :param string: str, string to write, each line is indented
        :param level: int, indentation level (>= 0)
        :param tag: str, optional XML tag
        :param indent_len: int, number of spaced to indent
        :param indent_char: chr, character to use in indenting
        :return: int, number of lines written
        -----------------------------------------------------------------------------------------"""
        nline = 0
        space = indent_char * level * indent_len
        if tag:
            fp.write("{}<{}>\n".format(space, tag))
            nline += 1

        for line in string.split('/n'):
            fp.write("{}{}\n".format(space, line))
            nline += 1

        if tag:
            fp.write("{}</{}>\n".format(space, tag))
            nline += 1

        return nline

    def XIOSwrite(self, fp):
        """-----------------------------------------------------------------------------------------
        Write the topology in XIOS XML format. The file comprises four sections
            <information>       metadata
            <stem_list>         a list of stems
            <edge_list>         list of edges (not needed)
            <adjacency>         adjacency matrix

        :param fp: file, a writeable file, e.g., sys.stdout
        :return:
        -----------------------------------------------------------------------------------------"""
        version = 2.1
        indent = 0
        indent_len = 2

        self.iwrite(fp, "<XIOS version='{}'>".format(version), 0)
        indent += indent_len
        space = ' ' * indent

        # information block
        self.iwrite(fp, '<information>', 1)

        if self.sequence:
            indent += indent_len
            space = ' ' * indent
            self.iwrite(fp, '<sequence>', 2)
            self.iwrite(fp, self.sequence_id, 3, 'id')
            self.iwrite(fp, self.sequence_doc, 3, 'doc')
            self.iwrite(fp, self.sequence, 3, 'unformatted')
            self.iwrite(fp, '</sequence>', 2)

            self.iwrite(fp, '<comment_list>', 2)
            for comment in self.comment:
                self.iwrite(fp, comment, 3, 'comment')
            self.iwrite(fp, '</comment_list>', 2)

        self.iwrite(fp, '</information>', 1)

        # stemlist block
        self.iwrite(fp, '<stem_list>', 1)
        self.iwrite(fp, self.format_stem_list(), 1)
        self.iwrite(fp, '</stem_list>', 1)

        # edge_list block
        self.iwrite(fp, '<edge_list>', 1)
        self.iwrite(fp, self.format_edge_list(), 1)
        self.iwrite(fp, '</edge_list>', 1)

        # adjacency matrix block
        self.iwrite(fp, '<adjacency>', 1)
        self.iwrite(fp, self.format_adjacency(), 1)
        self.iwrite(fp, '</adjacency>', 1)

        self.iwrite(fp, "</XIOS>".format(version), 0)

        return True

    def XIOSread(self, filename):
        """-----------------------------------------------------------------------------------------
        Read a XIOS topology file. The file comprises three sections, a list of stems ( read by
        stemlistRead), a list of edges (not needed), and an adjacency matrix (adjacencyRead).

        :return: int, nujmber of stems read
        -----------------------------------------------------------------------------------------"""
        fp = None
        try:
            fp = open(filename, 'r')
        except (OSError, IOError) as err:
            sys.stderr.write('Topology::XIOSread - error opening input file ({})\n'.
                             format(filename))
            sys.stderr.write('\t{}\n'.format(err))
            exit(1)

        x = etree.parse(fp)
        # print(etree.tostring(x))
        for section in x.xpath('//XIOS/*'):
            # print('section {}\n{}'.format(section.tag, section.text))
            if section.tag == 'information':
                self.parse_information(section)

            elif section.tag == 'stem_list':
                self.parse_stem_list(section.text)


            elif section.tag == 'edge_list':
                self.parse_edge_list(section.text)

            elif section.tag == 'adjacency':
                self.parse_adjacency(section.text)

            else:
                sys.stderr.write('Topology::XIOSread - unknown XML tag ({})\n'.format(section.tag))

        return 1

    def parse_information(self, x, clear=True):
        """-----------------------------------------------------------------------------------------
        Parse the information section of the XIOS formatted topology file.  Information  may
        include the sequence, provenence, or other metadata.  Information has the
        following structure, all elements are optional.

        All information fields have white space stripped from begin and end

        <information>
            <comment_list>
                <comment></comment>
            </comment_list>
            <sequence>
                <id></id>
                <doc></doc>
                <unformatted></unformatted>
            </sequence>

        
        :param x: str
        :param: clear: logical.  not implemented should clar the information block
        :return: True
        -----------------------------------------------------------------------------------------"""
        # if clear:
        #     info.clear()

        # field =  x.xpath('//information/*')
        # pass
        for element in x.iterfind('*'):
            if element.tag == 'comment_list':
                for comment in element.findall('*'):
                    if comment.tag == 'comment':
                        self.comment.append(comment.text.strip())
                        # print('\tcomment:\n{}'.format(comment.text))
                    else:
                        sys.stderr.write('Topology::parse_information - ')
                        sys.stderr.write('unknown comment_list element ({})\n'.format(element.tag))

            elif element.tag == 'sequence':
                for seqfield in element.findall('*'):
                    if seqfield.tag == 'id':
                        self.sequence_id = seqfield.text.strip()
                        # print('sequence id: {}'.format(seqfield.text))
                    elif seqfield.tag == 'doc':
                        self.sequence_doc = seqfield.text.strip()
                        # print('sequence doc: {}'.format(seqfield.text))
                    elif seqfield.tag == 'unformatted':
                        self.sequence = seqfield.text.strip()
                        # print('sequence text: {}'.format(seqfield.text))
                    else:
                        sys.stderr.write('Topology::parse_information - ')
                        sys.stderr.write('unknown sequence element ({})\n'.format(seqfield.tag))

            else:
                sys.stderr.write('Topology::parse_information - unknown information element ({})'.
                                 format(element.tag))

        return True

    def parse_stem_list(self, text, clear=True):
        """-----------------------------------------------------------------------------------------
        Parse the stem_list section of the XIOS formatted topology file. The stem list gives the
        location of each stem, and optionally, the viennna strings for each half stem. By default
        this method does not clear the current contents of the stem, use clear=False if you
        don't want to delete the current stemlist

        <stem_list>
             0 185.5 [   4  14   357 368 ]     (((((((((((   ).))))))))))
             1 176.5 [  15  21   332 338 ]         (((((((   )))))))
             2  35.0 [  23  32    38  48 ]      ((((((((((   ).)))))))))
             3 204.5 [  53  61   348 355 ]       (((.(((((   ))))))))
             4 148.5 [  62  65   232 235 ]            ((((   ))))
             5 170.0 [  69  72   268 271 ]            ((((   ))))
             6 152.5 [  74  78   227 231 ]           (((((   )))))
             7  86.0 [  79  83    89  93 ]           (((((   )))))
             8 100.5 [  94  98   103 108 ]           (((((   )))).)
             9 168.0 [ 109 118   218 226 ]      ((...(((((   )))..))))
            10 146.0 [ 131 144   148 165 ]  ((((((((((((((   )))))).))).)))..))
            11 182.0 [ 173 178   186 191 ]          ((((((   ))))))
            12 202.5 [ 192 200   205 215 ]       (((((((((   ))).))))).)
            13 267.0 [ 239 242   292 295 ]            ((((   ))))
            14 267.5 [ 249 254   281 286 ]          ((((((   ))))))
            15 269.0 [ 262 266   272 276 ]           (((((   )))))
            16 311.5 [ 302 309   314 321 ]        ((((((((   ))))))))
        </stem_list>

        description
                16 311.5 [ 302 309   314 321 ]        ((((((((   ))))))))
        field   0  1     2 3   4     5   6   7        8          9
        0   stem name                       str
        1   stem center                     float
        2   delimiter                       chr
        3   left half stem begin            int
        4   left half stem end              int
        5   right half stem begin           int
        6   right half stem end             int
        7   delimiter                       chr
        8   left half stem Vienna string    str optional
        9   right half stem Vienna string   str optional

        :param text: str
        :return: True
        -----------------------------------------------------------------------------------------"""
        stems = self.stem_list
        if clear:
            stems.clear()

        for line in text.split('\n'):

            if not line.strip():
                continue

            stem = {}
            field = line.split()
            stem['name'] = field[0]
            stem['center'] = float(field[1])
            stem['left_begin'] = int(field[3])
            stem['left_end'] = int(field[4])
            stem['right_begin'] = int(field[5])
            stem['right_end'] = int(field[6])
            if len(field) > 7:
                stem['left_vienna'] = field[8]
                stem['right_vienna'] = field[9]
            stems.append(stem)

        return True

    def parse_edge_list(self, text):
        """-----------------------------------------------------------------------------------------
        The edge_list section is for human reading purposes only.  I contains the same
        information, and is written from, the adjacency matrix.

        <edge_list>
             0:  1i  2i  3i  4i  5i  6i  7i  8i  9i 10i 11i 12i 13i 14i 15i 16i
             1:  2i  3o  4i  5i  6i  7i  8i  9i 10i 11i 12i 13i 14i 15i 16i
             2:
             3:  4i  5i  6i  7i  8i  9i 10i 11i 12i 13i 14i 15i 16i
             4:  5o  6i  7i  8i  9i 10i 11i 12i
             5:  6i  7i  8i  9i 10i 11i 12i 13o 14o 15o
             6:  7i  8i  9i 10i 11i 12i
             7:
             8:
             9: 10i 11i 12i
            10:
            11:
            12:
            13: 14i 15i
            14: 15i
            15:
            16:
        </edge_list>

        :param text: str
        :return: True
        -----------------------------------------------------------------------------------------"""

        return True

    def parse_adjacency(self, text):
        """-----------------------------------------------------------------------------------------
        Parse the adjacency section of the XIOS formatted topology file. If the section is absent,
        return False, otherwise True

        :param text: str
        :return: logical
        -----------------------------------------------------------------------------------------"""
        n = 0
        cols = 0
        for line in text.split('\n'):

            if not line.strip():
                continue

            n += 1
            if n == 1:
                # header line
                cols = int(line.split()[-1]) + 1
                adjacency = [['-' for _ in range(cols)] for _ in range(cols)]
                continue

            # line of adjacency matrix
            element = line.split()
            for i in range(1, cols + 1):
                adjacency[n - 2][i - 1] = element[i]

        self.adjacency = adjacency
        return True

    @staticmethod
    def sample(self, n):
        """-----------------------------------------------------------------------------------------
        randomly sample a connected subgraph of size=n from the topology.  The sampled graph will be
        connected by non-s edges, but size may be less than n if the graph being sampled is smaller
        than size or has disconnected segmentsSampling is based on the adjacency matrix.

        :param n:
        :return: topology
        -----------------------------------------------------------------------------------------"""
        # randomly determine starting vertex
        adj = self.adjacency
        nvertex = len(adj)

        vlist = []
        neighbor = []
        v0 = random.randrange(nvertex)
        size = 0

        while size < n:

            vlist.append(v0)
            size += 1

            # update list of neighbors
            for v1 in range(nvertex):
                if adj[v0][v1] == 's' or v1 == v0:
                    # skip s edges and self
                    continue

                if v1 in vlist:
                    # don't add already selected vertices to neighbors
                    continue

                if v1 not in neighbor:
                    # existing neighbor
                    neighbor.append(v1)

            if len(neighbor) == 0:
                # if there are no neighbors, you must stop
                break

            # select new vertex and remove from current neighbor list
            v0 = random.choice(neighbor)
            neighbor.remove(v0)

        # end of vertex selection loop
        vlist.sort()

        return vlist

    @staticmethod
    def samplebyweight(self, n, w):
        """-----------------------------------------------------------------------------------------
        randomly sample a connected subgraph of size=n from the topology weighted by w.  The sampled
        graph will be
        connected by non-s edges, but size may be less than n if the graph being sampled is smaller
        than size or has disconnected segmentsSampling is based on the adjacency matrix.

        :param n:
        :return: topology
        -----------------------------------------------------------------------------------------"""
        # randomly determine starting vertex
        adj = self.adjacency
        nvertex = len(adj)

        vlist = []
        neighbor = []
        idx = [i for i in range(nvertex)]
        v0 = random.choices(idx, weights=w)[0]
        size = 0

        while size < n:

            vlist.append(v0)
            size += 1

            # update list of neighbors
            for v1 in range(nvertex):
                if adj[v0][v1] == 's' or v1 == v0:
                    # skip s edges and self
                    continue

                if v1 in vlist:
                    # don't add already selected vertices to neighbors
                    continue

                if v1 not in neighbor:
                    # existing neighbor
                    neighbor.append(v1)

            if len(neighbor) == 0:
                # if there are no neighbors, you must stop
                break

            # select new vertex and remove from current neighbor list
            ww = [w[i] for i in neighbor]
            v0 = random.choices(neighbor,weights=ww)[0]
            neighbor.remove(v0)

        # end of vertex selection loop
        vlist.sort()

        return vlist

    def sample_topology(self, n):
        """-----------------------------------------------------------------------------------------
        Return a new topology sampled from the current one.
        :param n: int, number of stems to sample
        :param sample:
        :return: Topology
        -----------------------------------------------------------------------------------------"""
        vlist = self.sample(self, n)

        # build the new topology from the current one
        newtopo = Topology()
        for s in self.stem_list:
            if int(s['name']) in vlist:
                newtopo.stem_list.append(s)

        adj = self.adjacency
        newadj = newtopo.adjacency
        for row in vlist:
            newadj.append([])
            for v in vlist:
                newadj[-1].append(adj[row][v])

        # copy the old sequence information, add a new tracking comment
        newtopo.sequence = self.sequence
        newtopo.sequence_id = self.sequence_id
        newtopo.sequence_doc = self.sequence_doc
        newtopo.comment = copy.copy(self.comment)
        now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        newtopo.comment.append('history {} creator topology::sample(n={})'.format(now, n))

        return newtopo

    def sample_xios(self, n):
        """-----------------------------------------------------------------------------------------
        Return a xios structure sampled from the current topology.

        :param n: int, number of stems to sample
        :param sample:
        :return: Xios object (see xios.py)
        -----------------------------------------------------------------------------------------"""
        edge = {'i':0, 'j':1, 'o':2, 's':3, 'x':4 }
        vlist = self.sample(self, n)

        adj = self.adjacency
        struct = []
        for row in vlist:
            for col in vlist:
                if col <= row:
                    continue
                if adj[row][col] in ('i', 'j', 'o'):
                    struct.append([row, col, edge[adj[row][col]]])

        return Xios(list=struct)

    def sample_xios_weighted(self, n, w):
        """-----------------------------------------------------------------------------------------
        Return a xios structure sampled from the current topology.

        :param n: int, number of stems to sample
        :param w: list of weights
        :param sample:
        :return: Xios object (see xios.py)
        -----------------------------------------------------------------------------------------"""
        edge = {'i':0, 'j':1, 'o':2, 's':3, 'x':4 }
        vlist = self.samplebyweight(self, n, w)

        adj = self.adjacency
        struct = []
        for row in vlist:
            for col in vlist:
                if col <= row:
                    continue
                if adj[row][col] in ('i', 'j', 'o'):
                    struct.append([row, col, edge[adj[row][col]]])

        return Xios(list=struct), vlist


####################################################################################################
####################################################################################################
class SerialRNA(list):
    """=============================================================================================
    for working with RNAS encoded for example as 001212, meaning ( ) ( [ ) ]
    ============================================================================================="""

    def connected(self):
        """-----------------------------------------------------------------------------------------
        return the connected graph(s) in the current structure. The entire graph is connected if one
        structure is returned

        :return: list of SerialRNA
        -----------------------------------------------------------------------------------------"""
        component = []
        open = []
        begin = 0
        for pos in range(len(self)):
            if self[pos] in open:
                open.remove(self[pos])
                if len(open) == 0:
                    component.append(SerialRNA(self[begin:pos + 1]))
                    begin = pos + 1
            else:
                open.append(self[pos])

        if len(component) == 1:
            return [self]
        else:
            return component

    def addstemall(self):
        """-----------------------------------------------------------------------------------------
        return the set of structures with one additional stem explicitly added at all possible
        positions. This is just to check the other enumeration approaches.

        :return: list of SerialRNA
        -----------------------------------------------------------------------------------------"""
        # make a new structure with the stem number incremented by 1
        newlen = len(self) + 2
        children = []

        for begin in range(0, newlen - 1):
            for end in range(begin + 1, newlen):
                extended_rna = [None for _ in range(newlen)]
                extended_rna[begin] = 0
                extended_rna[end] = 0
                newpos = 0
                for pos in self:
                    while extended_rna[newpos] is not None:
                        newpos += 1
                    extended_rna[newpos] = pos + 1

                children.append(SerialRNA(extended_rna))

        return children

    def addstemleft(self):
        """-----------------------------------------------------------------------------------------
        return the set of structures with one additional stem added at all possible position

        :return: list of SerialRNA
        -----------------------------------------------------------------------------------------"""
        # make a new structure with the stem number incremented by 1
        # base = []
        newlen = len(self) + 2
        half = newlen // 2

        children = []
        # for pos in self:
        #     base.append(pos + 1)

        for begin in range(0, half - 1):
            for end in range(begin + 1, newlen):
                extended_rna = [None for _ in range(len(self) + 2)]
                extended_rna[begin] = 0
                extended_rna[end] = 0
                newpos = 0
                for pos in self:
                    while extended_rna[newpos] is not None:
                        newpos += 1
                    extended_rna[newpos] = pos + 1

                children.append(SerialRNA(extended_rna))

        return children

    def addstemzero(self):
        """-----------------------------------------------------------------------------------------
        return the set of structures with one additional stem such that the beginning of the stem
        is at position zero

        :return: list of SerialRNA
        -----------------------------------------------------------------------------------------"""
        # make a new structure with the stem number incremented by 1
        base = []
        children = []
        newlen = len(self) + 2
        for pos in self:
            base.append(pos + 1)

        for end in range(1, newlen):
            extended_rna = [None for _ in range(newlen)]
            extended_rna[0] = 0
            extended_rna[end] = 0
            basepos = 0
            for pos in range(1, newlen):
                if extended_rna[pos] is None:
                    extended_rna[pos] = base[basepos]
                    basepos += 1

            children.append(SerialRNA(extended_rna))

        return children

    def subtractstem(self):
        """-----------------------------------------------------------------------------------------
        return a list of the graphs made by subtracting one stem.  Returned graphs are in canonical
        form and unique

        :return: list of SerialRNA
        -----------------------------------------------------------------------------------------"""
        unique = {}
        for stem_num in range(len(self) // 2):
            parent = SerialRNA(self)
            parent.remove(stem_num)
            parent.remove(stem_num)
            for connected in parent.connected():
                connected.canonical()
                if connected == [0, 0]:
                    continue
                unique[connected.tostring()] = connected

        return list(unique.values())

    def canonical(self):
        """-----------------------------------------------------------------------------------------
        convert graph to canonical form.  In canonical form the stems occur in increasing numerical
        order beginning at zero

        :return: True
        -----------------------------------------------------------------------------------------"""
        stem = 0
        map = {}
        for pos in self:
            if pos not in map:
                map[pos] = stem
                stem += 1

        for pos in range(len(self)):
            self[pos] = map[self[pos]]

        return True

    def reverse(self):
        """-----------------------------------------------------------------------------------------
        reverse the rna structure, numbering is not converted to canonical
        :return: True
        -----------------------------------------------------------------------------------------"""
        begin = 0
        end = len(self) - 1
        while begin < end:
            self[begin], self[end] = self[end], self[begin]
            begin += 1
            end -= 1

        return True

    def canonical_fbstr(self):
        """-----------------------------------------------------------------------------------------
        return the graph and its reverse as strings in canonical form

        :return: str, str, forward and backward canonical strings
        -----------------------------------------------------------------------------------------"""
        self.canonical()
        forward = self.tostring()

        self.reverse()

        self.canonical()
        backward = self.tostring()

        return forward, backward

    def tostring(self):
        """-----------------------------------------------------------------------------------------
        return a string representing the stucture.  the string is the concatenation of the digits.

        :return: string
        -----------------------------------------------------------------------------------------"""
        return ''.join(str(x) for x in self)

    def fromstring(self, string):
        """-----------------------------------------------------------------------------------------
        convert a string of digits, e.g., 010122, to a SerialRNA

        :param string:
        :return: int, length of structure
        -----------------------------------------------------------------------------------------"""
        for c in list(string):
            self.append(int(c))

        return len(self)


####################################################################################################
####################################################################################################
class PairRNA:
    """=============================================================================================

    Synopsis
        from graph import PairRNA

        graph = PairRNA()
        graph.fromList([0,1,1,0])
            or
        graph = PairRNA(l=[0,1,1,0])

        graph.reverse()     # reverse the order of stems left-right
    ============================================================================================="""

    def __init__(self, inlist=[]):
        """-----------------------------------------------------------------------------------------
        Internal data structure is a list of lists:
        each element in the list of the begin, end position of a stem.  The stems
        are the coordinates of the left and right half stems as integers along a line.

        two nested stems are [ [0,3], [1,2] ]
        a simple pseudoknot is [ [0,2], [1,3] ]
        -----------------------------------------------------------------------------------------"""
        self.pairs = []
        self.nstem = 0

        if inlist:
            self.from_SerialRNA(inlist)

    def __str__(self, sep='_'):
        """-----------------------------------------------------------------------------------------
        Return a serialized version of the pair structure
        :return: string
        -----------------------------------------------------------------------------------------"""

        s = ''
        for p in self.pairs:
            try:
                s += '({},{}) '.format(p[0], p[1])
            except IndexError:
                # undefined stem
                s += '(.,.) '

        return s.rstrip('_')

    def __len__(self):
        """-----------------------------------------------------------------------------------------
        length is the number of stems
        -----------------------------------------------------------------------------------------"""
        return len(self.pairs)

    def from_SerialRNA(self, g):
        """"----------------------------------------------------------------------------------------
        read a graph in list format as a list.  returns a list of lists with the begin/end position
        of each stem (pair format)
        :param g: list, structure in list format
        :return: integer, number of stems
        -----------------------------------------------------------------------------------------"""
        self.nstem = int(len(g) / 2)
        self.pairs = [[] for _ in range(self.nstem)]

        # the values indicate the stem number, the index indicates the position
        for i in range(len(g)):
            self.pairs[g[i]].append(i)

        return self.nstem

    def from_SerialRNA_string(self, g, sep=' '):
        """-----------------------------------------------------------------------------------------
        the input list is a string separated by sep, e.g. '0 1 1 0'

        :param g: string, input graph as a SerialRNA string
        :param sep: string, separation character in input string
        :return: integer, number of stems
        -----------------------------------------------------------------------------------------"""
        pairs = self.pairs
        pairs.clear()
        nstem = 0

        value = g.split(sep)
        for n in range(len(value)):
            i = int(value[n])
            while i >= nstem:
                # create stems if seen for the first time
                pairs.append([])
                nstem += 1
            pairs[i].append(n)

        return len(pairs)

    def to_SerialRNA(self):
        """-----------------------------------------------------------------------------------------
        return the structure in list format (s a list)
        :return g: list, structure in list format
        TODO convert to SerialRNA object
        -----------------------------------------------------------------------------------------"""
        g = [0 for _ in range(self.nstem * 2)]

        stem = 0
        for pair in self.pairs:
            g[pair[0]] = stem
            g[pair[1]] = stem
            stem += 1

        return g

    def from_vienna(self):
        """-----------------------------------------------------------------------------------------
        Read structure in Vienna format
        :return: integer, number of stems
        TODO implement
        -----------------------------------------------------------------------------------------"""
        pass

    def to_vienna(self, pad=' '):
        """-----------------------------------------------------------------------------------------
        return a string with the structure in vienna format.  This is basically the abstract shapes
        format with support for pseudoknots.

        :return: string
        TODO test
        -----------------------------------------------------------------------------------------"""
        bracket = [['(', ')'], ['[', ']'], ['{', '}'], ['<', '>'], [':', ':']]
        vienna = ['.' for _ in range(self.nstem * 2)]
        list_format = self.to_SerialRNA()

        level = []
        knot = True
        for stemi in self.pairs:
            for l in range(len((level))):
                knot = True
                for stem_num in range(len(level[l])):
                    stemj = self.pairs[stem_num]
                    if stemi[0] > stemj[1] or (stemi[0] > stemj[0] and stemi[1] < stemj[1]):
                        # stemi is nested or serial with the stem using this level (stemj)
                        knot = False
                    else:
                        knot = True
                        break

                if not knot:
                    level[l].append(stemi)
                    break
            if knot:
                level.append([])
                level[-1].append(stemi)

        for l in range(len(level)):
            for stem in level[l]:
                vienna[stem[0]] = bracket[l][0]
                vienna[stem[1]] = bracket[l][1]

        return pad.join(vienna)

    def reverse(self):
        """-----------------------------------------------------------------------------------------
        reverse the positions of the stems by converting to maxpos-pos and resorting in order
        of first coordinate

        :return: None
        TODO test
        -----------------------------------------------------------------------------------------"""
        m = self.nstem * 2 - 1
        for i in range(self.nstem):
            self.pairs[i][0], self.pairs[i][1] = m - self.pairs[i][1], m - self.pairs[i][0]

        self.pairs.sort(key=lambda k: k[0])

        return self.pairs

    def connected(self):
        """-----------------------------------------------------------------------------------------
        Returns True if graph is i-o connected
        :return: True / False

        TODO test
        -----------------------------------------------------------------------------------------"""
        pos = 1
        for d in self.depth():
            if pos == 0:
                continue
            pos += 1
            if d == 0:
                break

        if pos < len(self) * 2:
            return False

        return True

    def depth(self):
        """-----------------------------------------------------------------------------------------
        Return  list with the nesting depth at each position of the graph in list format.  This is
        useful for mountain plots and determining connectivity.  If the level reaches zero before
        the last position, the graph is disconnected.

        :return: list
        TODO test
        -----------------------------------------------------------------------------------------"""
        depth = []
        d = 0
        stem = []

        for s in self.to_SerialRNA():
            if s in stem:
                # stem seen before
                d -= 1
            else:
                # new stem
                stem.append(s)
                d += 1

            depth.append(d)

        return depth


####################################################################################################
####################################################################################################
class RNAstructure:
    """=============================================================================================
    RNA structure class is for working with output from the RNAStructure package
    a single RNA structure

    TODO needs testing after refactoring

    ============================================================================================="""

    def __init__(self):
        self.sequence = ''
        self.length = 0
        self.energy = None  # not defined for probknot
        self.pair = []  # base number of the paired base
        self.stemlist = []
        self.adjacency = None
        self.id = None
        self.filename = None

    def CTRead(self, filename):
        """-----------------------------------------------------------------------------------------
        Read in RNAstructure CT file

        format example:
          240  ENERGY = -49.3   mr_s129
            1 A       0    2    0    1
            2 A       1    3    0    2
            3 C       2    4    0    3
            4 C       3    5    0    4
            5 A       4    6    0    5

            base_number base previous_base next_base paired_base base_number

        usage
            rna.CTRead(filename)
        :param filename: string, filename with CT formatted structure
        :return: integer, number of bases
        -----------------------------------------------------------------------------------------"""
        nbase = 0
        path = filename.split('/')
        self.filename = path[-1]
        with open(filename, 'r') as ct:
            line = ct.readline()
            # print('firstline:', line)
            # print('field:',field)

            if line.find('ENERGY') >= 0:
                field = line.split()
                self.energy = field[3]

            else:
                # probknot file
                field = line.split()
                self.length = int(field[0])
                self.id = field[1]

            self.pair = [0] * (self.length + 1)

            for line in ct:
                n, base, prev, following, pair, n2 = line.split()
                # print('n:', n, 'base:', base, 'pref:', prev, 'next:', next, 'pair:', pair, 'n2:',
                #       n2)
                self.sequence += base
                if pair != '0':
                    self.pair[int(pair)] = int(n)
                    self.pair[int(n)] = int(pair)
                nbase += 1

        return nbase

    def __str__(self):
        """-----------------------------------------------------------------------------------------
        string representation of a structure
        :return: string
        -----------------------------------------------------------------------------------------"""
        rnastr = 'RNA:\n'
        for key in self.__dict__:
            rnastr += '{0} = {1}\n'.format(key, self.__dict__[key])

        return rnastr

    def stemlist_from_pairs(self, unpaired=2):
        """-----------------------------------------------------------------------------------------
        Construct the stemlist from the paired base list in stem.pair
        :return: integer, number of stems in stemlist
        -----------------------------------------------------------------------------------------"""
        maxgap = unpaired + 1
        nstem = 0
        instem = False
        for pos in range(0, len(self.pair) - 1):
            if self.pair[pos] == 0 or self.pair[pos] < pos:
                continue

            if instem:
                # currently in a stem
                lgap = pos - stem.lend - 1
                rgap = stem.rbegin - self.pair[pos] - 1
                if lgap >= maxgap or rgap >= maxgap:
                    # gap is too big, end old stem
                    stem.trimVienna()
                    instem = False
                else:
                    # extend current stem
                    stem.lend = pos
                    stem.rbegin = self.pair[pos]
                    stem.lvienna += '{}('.format('.' * lgap)
                    stem.rvienna = '){}'.format('.' * rgap) + stem.rvienna
                    continue

            # not in a stem, start a new stem
            stem = Stem()
            nstem += 1
            self.stemlist.append(stem)
            instem = True

            stem.lbegin = pos
            stem.rend = self.pair[pos]
            stem.lend = pos
            stem.rbegin = self.pair[pos]
            stem.lvienna = '('
            stem.rvienna = ')'

        return nstem

    def stemListGetFromGraph(self, g):
        """-----------------------------------------------------------------------------------------
        Create stemlist from an PairRNA object (graph.py).  Since the RNA graph is abstract the
        start and end locations of the stems correspond to the stem position
        :param g: PairRNA object
        :return: integer, number of stems in stemlist
        -----------------------------------------------------------------------------------------"""
        vienna = g.to_vienna().split()
        nstem = 0
        for stem in g.pairs:
            s = Stem()
            s.lbegin = stem[0]
            s.lend = stem[0]
            s.rbegin = stem[1]
            s.rend = stem[1]
            s.lvienna = vienna[stem[0]]
            s.rvienna = vienna[stem[1]]
            self.stemlist.append(s)
            nstem += 1

        self.nstem = nstem
        return nstem

    def stemlist_format(self):
        """-----------------------------------------------------------------------------------------
        Returns a string with the stemlist formatted according to Stem.formatted()
        :return: string
        -----------------------------------------------------------------------------------------"""
        n = 0
        stemstr = ''
        for stem in self.stemlist:
            n += 1
            stemstr += '{0}\t{1}\n'.format(n, stem.formatted())
        return stemstr

    def adjacencyGet(self):
        """-----------------------------------------------------------------------------------------
        Calculate an adjacency matrix from a stemlist. assumes stems are ordered by the beginning
        of the left half-stem.
        :return: dict, keys are edge types, values are counts
        -----------------------------------------------------------------------------------------"""
        nstem = len(self.stemlist)

        edges = {'i': 0, 'j': 0, 'o': 0, 's': 0, 'x': 0}
        a = [[0 for _ in range(nstem)] for _ in range(nstem)]

        for i in range(nstem):
            stem_i = self.stemlist[i]

            for j in range(i + 1, nstem):
                stem_j = self.stemlist[j]

                if stem_i.rend < stem_j.lbegin:
                    # serial edge
                    a[i][j] = 's'
                    a[j][i] = 's'
                    edges['s'] += 1

                elif stem_i.lend < stem_j.lbegin and stem_i.rend < stem_j.rbegin:
                    # overlap edge (pseudoknot)
                    a[i][j] = 'o'
                    a[j][i] = 'o'
                    edges['o'] += 1

                elif stem_i.lend < stem_j.lbegin and stem_i.rbegin > stem_j.rend:
                    # included edge (directed, j is inside i)
                    a[i][j] = 'i'
                    a[j][i] = 'j'
                    edges['i'] += 1

                else:
                    # excluded edge (stems overlap)
                    a[i][j] = 'x'
                    a[j][i] = 'x'
                    edges['x'] += 1

        self.adjacency = a
        return edges

    def adjacency_format(self):
        """-----------------------------------------------------------------------------------------

        :return: string, formatted version of adjacency matrix
        -----------------------------------------------------------------------------------------"""
        adjstr = ''
        # TODO set width based on matrix size?
        width = 4
        if not self.adjacency:
            adjstr = 'None'
            return adjstr

        size = len(self.adjacency)
        for i in range(size):
            for j in range(size):
                # must cast to str so zeros format correctly
                adjstr += '{:{}}'.format(str(self.adjacency[i][j]), width)
            adjstr += '\n'

        return adjstr

    def edgelist(self, include="ijo", whole=False):
        """-----------------------------------------------------------------------------------------
        An edgelist is an array of lists.  each row corresponds to a stem (vertex).  The values are
        tuples with the number and type of nodes with edges
        :return: list, edgelist
        -----------------------------------------------------------------------------------------"""
        elist = []
        if not self.adjacency:
            return elist

        size = len(self.adjacency)
        a = self.adjacency
        for i in range(size):
            e = []
            elist.append(e)
            begin = i + 1
            if whole:
                begin = 0

            for j in range(begin, size):
                if i == j:
                    continue
                if a[i][j] in include:
                    e.append([j, a[i][j]])
        return elist

    def edgelist_format(self, include='ijo', whole=False):
        """-----------------------------------------------------------------------------------------
        Return a formatted version of the edgelist
        :param include: string, list of edgetypes to include
        :param whole: boolean, If true return the square matrix, otherwise triangular
        :return: string
        -----------------------------------------------------------------------------------------"""
        edgestr = ''
        # e = self.edgelist(include,whole)
        n = 0
        for edge in self.edgelist(include, whole):
            n += 1
            edgestr += '{}: '.format(n)
            for neighbor in edge:
                edgestr += '{}{}, '.format(neighbor[0], neighbor[1])

            edgestr = edgestr.rstrip(', ')
            edgestr += '\n'

        return edgestr


####################################################################################################
####################################################################################################
class Stem:
    """=============================================================================================
    coordinates of stem regions (including allowed gaps) and the corresponding Vienna strings.
    ============================================================================================="""

    def __init__(self):
        """-----------------------------------------------------------------------------------------
        Stem constructor
        begin-end coordinates are small->big for left side
        end-begin for right side
        -----------------------------------------------------------------------------------------"""
        self.lbegin = 0
        self.lend = 0
        self.rbegin = 0
        self.rend = 0
        self.lvienna = ''
        self.rvienna = ''

    def formatted(self):
        """-----------------------------------------------------------------------------------------
        Return a formatted version of the stem:
            left begin pos
            left end pos
            right begin pos
            right end pos
            left Vienna string
            right Vienna string
        :return: string, with stem coordinates and Vienna structure
        -----------------------------------------------------------------------------------------"""
        return '{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format(self.lbegin, self.lend, self.rbegin, self.rend,
                                                     self.lvienna, self.rvienna)

    def trimVienna(self):
        """-----------------------------------------------------------------------------------------
        Removes leading and trailing unpaired positions (.) from vienna strings
        :return: True
        -----------------------------------------------------------------------------------------"""
        self.lvienna = self.lvienna.rstrip('.')
        self.rvienna = self.rvienna.lstrip('.')
        return True


####################################################################################################
# --------------------------------------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------------------------------------
####################################################################################################

if __name__ == '__main__':
    def test_pair():
        """-----------------------------------------------------------------------------------------
        PairRNA object test functions

        :return:
        -----------------------------------------------------------------------------------------"""

        # PairRNA
        print('Testing PairRNA')
        g = [0, 1, 1, 0]
        pair = PairRNA()
        pair.from_SerialRNA(g)
        print('Serial {} => {}'.format(g, pair))

        gstring = '0 1 1 0'
        pair.from_SerialRNA_string(gstring)
        print('Serial string {} => {}'.format(gstring, pair))

        gstring = '0,2,2,0'
        pair.from_SerialRNA_string(gstring, sep=',')
        print('Serial string {} => {}'.format(gstring, pair))

        return True


    def test_topology():
        """-----------------------------------------------------------------------------------------
        Topology object test functions
        -----------------------------------------------------------------------------------------"""
        print('\nTesting Topology')

        top = Topology()
        top.XIOSread('data/rnasep_a1.Buchnera_APS.xios')
        top.XIOSwrite(sys.stdout)

        # sample a subgraph and confirm it is independent
        sample = top.sample_topology(5)
        top.stem_list[0] = {}
        top.adjacency = []
        top.comment = []
        sample.XIOSwrite(sys.stdout)

        # sample xios from topology
        top2 = Topology()
        top2.XIOSread('data/rnasep_a1.Buchnera_APS.xios')
        x = top2.sample_xios(5)
        print(x)

        return True


    test_pair()
    test_topology()

    exit(0)

    # for rna in graphs:
    #     g = RNAGraph(rna)
    #     three = RNAstructure()
    #     three.stemListGetFromGraph(g)
    #     print('Stemlist\n')
    #     print(three.stemlist_format())
    #
    #     edges = three.adjacencyGet()
    #     print('\nAdjacency matrix\n')
    #     print(three.adjacency_format())
    #
    #     print('\nEdgelist\n')
    #     e = three.edgelist()
    #     print(three.edgelist_format())
    #
    # topology
    #
    #     # from PairRNA
    #
    #     print('\nTesting connectivity')
    #     graph = RNAGraph(inlist=[0, 0, 1, 1, 2, 2])
    #     print('    ', graph.pairs)
    #     if not graph.connected():
    #         print('Not connected')
    #
    #     print('\nlist format')
    #     structure = [0, 1, 2, 1, 0, 2]
    #     print('    input list', structure)
    #     graph = RNAGraph(inlist=structure)
    #     print('    serialized', str(graph))
    #     print('    pairs', graph.pairs)
    #     print('    list', graph.to_SerialRNA())
    #     print('    vienna', graph.to_vienna())
    #
    #     print('\n    reverse')
    #     print('    reversed', graph.reverse())
    #     print('    pairs', graph.pairs)
    #     print('    list', graph.to_SerialRNA())
    #     print('    vienna', graph.to_vienna())
    #
    #     print('\nenumerating: size, len, total')
    #     total = 0
    #     for size in range(1, 8):
    #         g = enumerateRNATopology(size)
    #         total += len(g)
    #         print('    ', size, len(g), total)
    #
    #     print('\n3 stem graphs')
    #     graphs = enumerateRNATopology(3)
    #     for g in graphs:
    #         pgraph = RNAGraph(g)
    #         print('\ngraph:', g)
    #         print('    pairs:', pgraph.pairs, end='\t=>\t')
    #         print('    list:', pgraph.to_SerialRNA())
    #         print('    vienna', pgraph.to_vienna())
    #
    #         pgraph.reverse()
    #         print('    reversed:', str(pgraph))
    #         print('    pairs:', pgraph.pairs, end='\t=>\t')
    #         glist = pgraph.to_SerialRNA()
    #         print('    reversed list:', glist)
    #         print('    vienna', pgraph.to_vienna())
    #
    #     # RNAstructure
    #
    #     rna = RNAstructure()
    #     rna.CTRead('data/mr_s129.probknot.ct')
    #     rna.stemlist_from_pairs()
    #     # print(rna)
    #     print('Stemlist\n')
    #     print(rna.stemlist_format())
    #
    #     edges = rna.adjacencyGet()
    #     print('edges', edges)
    #     print('\nAdjacency matrix\n')
    #     print(rna.adjacency_format())
    #
    #     print('\nEdgelist\n')
    #     e = rna.edgelist()
    #     for i in range(len(e)):
    #         print('{}:\t{}'.format(i, e[i]))
    #
    #     """ correct structure of mr_s129 by manual inspection
    #     29-31:231-229
    #     33-37:225-221
    #     41-44:220-217
    #     48-52:211-207
    #     55-58:68-65
    #     65-68:58-55
    #     76-80:192-188
    #     84-91:184-177
    #     95-98:155-152
    #     100-103:150-147
    #     106-108:145-143
    #     117-119:130-128
    #     129-130:119-117
    #     143-145:108-106
    #     147-150:103-100
    #     152-155:98-95
    #     163-165:176-174
    #     174-176:165-163
    #     177-184:91-84
    #     188-192:80-76
    #     201-205:216-212
    #     207-211:52-48
    #     212-216:205-201
    #     217-225:44-33
    #     229-231:31-29
    #     """

    exit(0)
