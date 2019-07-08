import sys
from lxml import etree


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

    many topologies will have only some of these
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
        print(etree.tostring(x))
        for section in x.xpath('//XIOS/*'):
            print('section {}\n{}'.format(section.tag, section.text))
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

        
        :param text: str
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
                        self.sequence_id = seqfield.text.strip()
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

        n = 0
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


class RNAstructure:
    """=============================================================================================
    RNA structure class is for working with output from the RNASTructure package
    a single RNA structure

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

    def stemListGetFromPairs(self, unpaired=2):
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
        Create stemlist from an RNAGraph object (graph.py).  Since the RNA graph is abstract the
        start and end locations of the stems correspond to the stem position
        :param g: RNAGraph object
        :return: integer, number of stems in stemlist
        -----------------------------------------------------------------------------------------"""
        vienna = g.toVienna().split()
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

    def stemlistFormat(self):
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

    def adjacencyFormat(self):
        """-----------------------------------------------------------------------------------------

        :return: string, formatted version of adjacency matrix
        -----------------------------------------------------------------------------------------"""
        adjstr = ''
        width = 4  # TODO set based on matrix size?
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

    def edgelistFormat(self, include='ijo', whole=False):
        """-----------------------------------------------------------------------------------------
        Return a formatted version of the edgelist
        :param include: string, list of edgetypes to include
        :param whole: boolean, If true return the square matrix, otherwise triangular
        :return: string
        -----------------------------------------------------------------------------------------"""
        edgestr = ''
        # e = self.edgelist(include,whole)
        n = 0
        for edge in self.edgelist(include,whole):
            n += 1
            edgestr += '{}: '.format(n)
            for neighbor in edge:
                edgestr += '{}{}, '.format(neighbor[0], neighbor[1])

            edgestr = edgestr.rstrip(', ')
            edgestr += '\n'

        return edgestr


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
            right begin pso
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


# ==================================================================================================
# main/test
# ==================================================================================================

if __name__ == '__main__':
    rna = RNAstructure()
    rna.CTRead('data/mr_s129.probknot.ct')
    rna.stemListGetFromPairs()
    # print(rna)
    print('Stemlist\n')
    print(rna.stemlistFormat())

    edges = rna.adjacencyGet()
    print('edges', edges)
    print('\nAdjacency matrix\n')
    print(rna.adjacencyFormat())

    print('\nEdgelist\n')
    e = rna.edgelist()
    for i in range(len(e)):
        print('{}:\t{}'.format(i, e[i]))


""" correct structure of mr_s129 by manual inspection
29-31:231-229
33-37:225-221
41-44:220-217
48-52:211-207
55-58:68-65
65-68:58-55
76-80:192-188
84-91:184-177
95-98:155-152
100-103:150-147
106-108:145-143
117-119:130-128
129-130:119-117
143-145:108-106
147-150:103-100
152-155:98-95
163-165:176-174
174-176:165-163
177-184:91-84
188-192:80-76
201-205:216-212
207-211:52-48
212-216:205-201
217-225:44-33
229-231:31-29
"""

from graph import enumerateRNATopology, RNAGraph
graphs = enumerateRNATopology(3)

for rna in graphs:
    g = RNAGraph(rna)
    three = RNAstructure()
    three.stemListGetFromGraph(g)
    print('Stemlist\n')
    print(three.stemlistFormat())

    edges = three.adjacencyGet()
    print('\nAdjacency matrix\n')
    print(three.adjacencyFormat())

    print('\nEdgelist\n')
    e = three.edgelist()
    print(three.edgelistFormat())

    top = Topology()
    top.XIOSread('data/rnasep_a1.Buchnera_APS.xios')

    exit(0)
