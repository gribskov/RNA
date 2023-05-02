"""=================================================================================================
    There are two main classes in this file.
    Stem object contain information about a single stem. The begin and end positions of the left
    and right halves of the stem and the Vienna (dot-bracket) structure.
    ================================================================================================
    Topology describes topologies are sets of stems that may include stems that are mutually
    exclusive (or even identical).RNAstructure, PairRNA, and SerialRNA describe single structures.

    Topology is the base class and provides functions for reading, parsing, and writing the XIOS
    format.  Topology is the ultimate arbiter of the XIOS format

    RNAstructure reads structures from CT files produced by the RNAstructure package (or by
    mfold/unafold).  Multiple structures can be read to create a single topology if the CT file is
    produced by the Fold program

    TODO: check if PairRNA and SerialRNA completely populate a usable Topology object
    PairRNA produces Topology objects from the pair representation (see below)
    SerialRNA produces Topology objects from the serial format (see below)

    ================================================================================================

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
                  1 8 2 3 4 5 6 7               (represented internall in serialized pair form)

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
        mfold and RNAStructure

    RNAml (currently not supported)
================================================================================================="""
import sys
import copy
from datetime import datetime
import random
from lxml import etree


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

    def __init__(self, *args, **kwds):
        """-----------------------------------------------------------------------------------------

        -----------------------------------------------------------------------------------------"""
        self.name = ''
        self.stem_list = []
        self.edge_list = []
        self.adjacency = []
        self.vienna = ''
        self.sequence = ''
        self.sequence_id = ''
        self.sequence_doc = ''
        self.sequence_length = 0
        self.comment = []

        for key in kwds:
            if key == 'xml' or key == 'xios':
                self.XIOSread(kwds[key])
            else:
                sys.stderr.write('Topology::init - unknown keyword ({})'.format(key))

    def format_stem_list(self):
        """-----------------------------------------------------------------------------------------
        Construct a formatted string in which the stemlist columns line up nicely. is a Stem
        object with attributes
            name
            center
            legin
            lend
            rbegin
            rend
            lvienna
            rvienna
        
        :return: str
        -----------------------------------------------------------------------------------------"""
        attrs = ['name', 'lbegin', 'lend', 'lvienna', 'rvienna', 'rbegin', 'rend']
        string = ''
        colmax = {}
        for s in self.stem_list:
            # for formatting, find the width of each attribute in the stem_list
            for column in attrs:
                colstr = '{}'.format(getattr(s, column))
                if column in colmax:
                    colmax[column] = max(colmax[column], len(colstr))
                else:
                    colmax[column] = len(colstr)

        fmt = ' {{:>{}}}  {{:{}}}  [ {{:{}}} {{:{}}} {{:{}}} {{:{}}} ]  {{:>{}}}  {{:{}}}\n'. \
            format(colmax['name'], colmax['rend'] + 2,
                   colmax['lbegin'], colmax['lend'],
                   colmax['rbegin'], colmax['rend'],
                   colmax['lvienna'], colmax['rvienna'])

        for s in self.stem_list:
            # stem center is average of lbegin and rend
            string += fmt.format(s.name, round((s.lbegin + s.rend) / 2, 1),
                                 s.lbegin, s.lend,
                                 s.rbegin, s.rend,
                                 s.lvienna, s.rvienna)

        return string.rstrip()

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
            string += '\n'
            r += 1

        return string.rstrip()

    def format_adjacency(self, fieldwidth=3):
        """-----------------------------------------------------------------------------------------
        Return a formatted string with the adjacency matrix.

        :param fieldwidth: int, with of field in formatted string
        :return: str
        -----------------------------------------------------------------------------------------"""
        string = ' ' * fieldwidth
        adj = self.adjacency

        size = max(len(str(len(adj))) + 1, fieldwidth)
        fmt = '{{:>{}}}'.format(size)

        string += ''.join([fmt.format(i) for i in range(len(adj))])
        r = 0
        for row in adj:
            string += '\n'
            string += fmt.format(r)
            string += ''.join([fmt.format(row[i]) for i in range(len(adj))])
            r += 1

        return string.rstrip()

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

        for line in string.split('\n'):
            if line:
                fp.write("{}{}\n".format(space, line))
            nline += 1

        if tag:
            fp.write("{}</{}>\n".format(space, tag))
            nline += 1

        return nline

    def edgelist_from_adjacency(self, include="ijo", whole=False):
        """-----------------------------------------------------------------------------------------
        An edgelist is an array of lists.  each row corresponds to a stem (vertex).  The values are
        tuples with the number and type of nodes with edges.  This function populates the edge_list
        of the topology object

        :param include: string, edge types to include
        :param whole: boolean, include self edges, why??
        :return: int, number of edges
        -----------------------------------------------------------------------------------------"""
        elist = []
        if not self.adjacency:
            return 0

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

        self.edge_list = elist

        return len(self.edge_list)

    def XIOSwrite(self, fp):
        """-----------------------------------------------------------------------------------------
        Write the topology in XIOS XML format. The file comprises four sections
            <information>       metadata
            <stem_list>         a list of stems
            <edge_list>         list of edges (not needed)
            <adjacency>         adjacency matrix

        :param fp: file, a writeable file, e.g., sys.stdout
        :return: True
        -----------------------------------------------------------------------------------------"""
        version = 2.1
        indent = 0
        indent_len = 2

        self.iwrite(fp, "<XIOS version='{}'>".format(version), 0)
        indent += indent_len
        # space = ' ' * indent

        # information block
        self.iwrite(fp, '<information>', 1)

        if self.sequence:
            indent += indent_len
            # space = ' ' * indent
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

    def XIOSread(self, file):
        """-----------------------------------------------------------------------------------------
        Read a XIOS topology file. The file comprises three sections, a list of stems ( read by
        stemlistRead), a list of edges (not needed), and an adjacency matrix (adjacencyRead).

        :return: int, number of stems read
        -----------------------------------------------------------------------------------------"""
        fp = None
        if isinstance(file, str):
            # file argument is string, try to open
            try:
                fp = open(file, 'r')
            except (OSError, IOError) as err:
                sys.stderr.write('Topology::XIOSread - error opening input file ({})\n'.
                                 format(file))
                sys.stderr.write('\t{}\n'.format(err))
                exit(1)
        else:
            # file is not str, assume it is a file pointer
            fp = file

        x = etree.parse(fp)
        # print(etree.tostring(x))
        for section in x.xpath('//XIOS/*'):
            # print('section {}\n{}'.format(section.tag, section.text))
            if section.tag == 'information':
                self.parse_information(section)

            elif section.tag == 'stem_list':
                self.parse_stem_list(section.text)

            elif section.tag == 'edge_list':
                pass
                # self.parse_edge_list(section.text)

            elif section.tag == 'adjacency':
                self.parse_adjacency(section.text)

            else:
                sys.stderr.write('Topology::XIOSread - unknown XML tag ({})\n'.format(section.tag))

        self.edgelist_from_adjacency()
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
        :param clear: logical.  not implemented should clear the information block
        :return: True
        -----------------------------------------------------------------------------------------"""
        if clear:
            pass
            # info.clear()

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
        :param clear, logical - clear the stems list
        :return: True
        -----------------------------------------------------------------------------------------"""
        stems = self.stem_list
        if clear:
            stems.clear()

        for line in text.split('\n'):

            if not line.strip():
                continue

            stem = Stem()
            field = line.split()
            stem.name = field[0]
            stem.center = float(field[1])
            stem.lbegin = int(field[3])
            stem.lend = int(field[4])
            stem.rbegin = int(field[5])
            stem.rend = int(field[6])
            if len(field) > 8:
                stem.lvienna = field[8]
                stem.rvienna = field[9]
            stems.append(stem)

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

    def stemlist_merge_case2(self, stem_center_threshold=2):
        """-----------------------------------------------------------------------------------------
          |0aaaaaaaa1|      |2aaaaaaaa3|
                |0bbbb1|   |2bbbb3|
    
        merges when
        1) the left half stem of the parent partially/completely overlaps the left
           halfstem of the child, AND
        2) the right half stem of the parent partially/completely overlaps left
           halfstem of the child, AND
        3) the stem center difference is less than stem_center_threshold
    
        in this situation, the child has some bases that lie outside the parent half
        stems on the right are merged by stemlist_deduplicate()
    
        a is the parent stem and b is the child, merging occurs when
        
        b.lbegin >= a.lbegin and b.lbegin <= a.lend and b.lend < a.rbegin
    
        ( stemcenter_diff <= stem_center_threshold )                            and
        ( b.lbegin >= a.lbegin and b.lbegin <= a.lend and b.lend < a.rbegin )   and
        ( b.rend >= a.rbegin and b.rend <= a.rend and b.rbegin>a.lend )         and
        ( b.lend > a.lend or b.rbegin < a.rbegin )                              and

        :param stem_center_threshold: int, allowed difference in center of merged stems 
        :return: 
        -----------------------------------------------------------------------------------------"""
        nmerged = 0
        for i in range(len(self.stem_list)):
            s1 = self.stem_list[i]
            for j in range(i + 1, len(self.stem_list)):
                s2 = self.stem_list[j]
                if s2.lbegin < s1.lend:
                    # stem 2 is to the right of stem 1, since stems are sorted by lbegin, no other
                    # stems can overlap
                    break

                if s2.rbegin > s1.rbegin or s2.rend > s1.rend or s2.rend < s1.rbegin:
                    # right half stem does not meet conditions of mergecase
                    break
                if abs(s1.center - s2.center) > stem_center_threshold:
                    break

                # qualifying overlap

    #     for my i ( nodelist ) {
    #     # skip already merged stems and the root node
    #     next if ( i->{id} =~ / x / | | i->{id} == 0 )
    #
    #     # the parent stem
    #     my ( a0, a1, a2, a3, ac ) = ( i->{pos}->[0], i->{pos}->[1],
    #     i->{pos}->[2], i->{pos}->[3],
    #     i->{pos}->[4]
    #     )
    #
    #     my children = {i->{children}}
    #     foreach my child ( children ) {
    #     next if ( child->{id} =~ / x / )  # skip already merged
    #     # the child stem
    #     my ( b0, b1, b2, b3, bc ) = ( child->{pos}->[0], child->{pos}->[1],
    #     child->{pos}->[2], child->{pos}->[3],
    #     child->{pos}->[4]
    #     )
    #
    #     if ( (b0 >= a0 & & b0 <= a1 & & b1 < a2) & &
    #     (b3 >= a2 & & b3 <= a3 & & b2 > a1) & &
    #     (b1 > a1 | | b2 < a2)                   ) {
    #     # condition for case 2, has been met,
    #     # Is there any need to test for unmatched bases? I think not...
    #     # TODO: check that the loop region is large enough to be a real loop?
    #
    #     my stemcenter_diff = abs(bc - ac)
    #     if ( stemcenter_diff <= stem_center_threshold ) {
    #     printf STDERR "   Case 2:Child node %s [ %d %d %d %d ] merged into ",
    #     child->{id}, b0, b1, b2, b3 unless quiet
    #
    #     printf STDERR "%s [ %d %d %d %d ]\n",
    #     i->{id}, a0, a1, a2, a3, ac  unless quiet
    #
    #     treeMergeNodes( i, child )
    #     nmerged++
    #     }
    #     }
    #     }
    #     }
    #     return nmerged
    #     }
    #
    #     # end of mergeCase2
    #
    #     # ------------------------------------------------------------------------------
    #     # mergeCase3
    #     #
    #     #   |0aaaaaaaa1|                    |2aaaaaaaa3|
    #     #                |0bbbb1|  |2bbbb3|
    #     #
    #     # merges child and parent nodes when
    #     # 1) the child is completely nested within the parent stem, AND
    #     # 2) the gap between child and parent stems is not greater than $gap_limit, AND
    #     # 3) the centers are less than $STEM_CENTER_THRESHOLD apart
    #     #
    #     # ( $b0 > $a1 && $b1 < $a2 && $b2 > $a1 && $b3 < $a2 ) and
    #     #                 $left_halfstem_gap <= $threshold and
    #     #                 $right_halfstem_gap <= $threshold and
    #     #                 $stemcenter_diff <= $STEM_CENTER_THRESHOLD
    #     #
    #     # ( $b0 > $a1 && $b1 < $a2 && $b2 > $a1 && $b3 < $a2 ) &&
    #     # ( $b0-$a1 <= $threshold && $a2-$b3 <= $threshold )   &&
    #     # ( $abs($bc - $ac) <= $STEM_CENTER_THRESHOLD )
    #     #
    #     # $b0 - $a1 is the gap between the parent and child left half stems
    #     # $a2 - $b3 is the gap between the parent and child right half stems
    #     # abs($bc - $ac) is the difference in the center position of the two stems
    #     #
    #     # USAGE
    #     #   $nmerged = mergeCase3( $gap_limit, $nodelist );
    #     #   Note: uses global variable $STEM_CENTER_THRESHOLD
    #     # ------------------------------------------------------------------------------
    #     sub
    #     mergeCase3
    #     {
    #     my( $threshold,
    #
    #     @nodelist
    #
    #     ) =
    #
    #     @_
    #
    #     ;
    #
    #     my $nmerged = 0;
    #     for my $i ( @ nodelist ) {
    #     # skip already merged stems and the root node
    #     next if ( $i->{id} =~ / x / | | $i->{id} == 0 );
    #
    #     # the parent stem
    #     my ( $a0, $a1, $a2, $a3, $ac ) = ( $i->{pos}->[0], $i->{pos}->[1],
    #     $i->{pos}->[2], $i->{pos}->[3],
    #     $i->{pos}->[4]
    #     );
    #
    #     my @ children = @ {$i->{children}};
    #     foreach my $child ( @ children ) {
    #     next if ( $child->{id} =~ / x / );
    #     my ( $b0, $b1, $b2, $b3, $bc ) = ( $child->{pos}->[0], $child->{pos}->[1],
    #     $child->{pos}->[2], $child->{pos}->[3],
    #     $child->{pos}->[4]
    #     );
    #
    #     if ( $b0 > $a1 & & $b1 < $a2 & & $b2 > $a1 & & $b3 < $a2 ) {
    #     # don't need to check that a loop is physically reasonable since
    #     # the nested stem has a loop
    #
    #     my $stemcenter_diff    = abs($bc - $ac);
    #     my $left_halfstem_gap  = $b0 - $a1;
    #     my $right_halfstem_gap = $a2 - $b3;
    #
    #     if ( $left_halfstem_gap <= $threshold and
    #     $right_halfstem_gap <= $threshold and
    #     $stemcenter_diff <= $STEM_CENTER_THRESHOLD ) {
    #
    #     printf STDERR "   Case 3: Child node %s [ %d %d %d %d ] merged into ",
    #     $child->{id}, $b0, $b1, $b2, $b3 unless $quiet;
    #
    #     printf STDERR "%s [ %d %d %d %d ]\n",
    #     $i->{id}, $a0, $a1, $a2, $a3, $ac  unless $quiet;
    #
    #     treeMergeNodes( $i, $child );
    #     $nmerged++;
    #     }
    #     }
    #     }
    #     }
    #     return $nmerged;
    #
    # }
    #
    # # end of mergeCase3
    #
    # # ------------------------------------------------------------------------------
    # # mergeCase4
    # #
    # #       |0aaaaaaaa1|           |2aaaaaaaa3|
    # #        |0bbbb1|     |2bbbb3|
    # #
    # #       |0aaaaaaaa1|           |2aaaaaaaa3|
    # #                    |0bbbb1|     |2bbbb3|
    # #
    # # merges nested child node with parent nodes when
    # # 1) one child halfstem overlaps partially/completely with the parent AND
    # # 2) the other half stem is nested inside the parent half stems AND
    # # 3) the nested half stem is less than $gap_limit away from the boundaries of
    # #    the parent AND
    # # 4) the stem center difference is less than $STEM_CENTER_THRESHOLD
    # #
    # # ( $b0 <= $a1 && $b0 >= $a0 && $b3 <  $a2 && $b2 >  $a1 ) ||
    # # ( $b0 > $a1  && $b1 < $a2  && $b3 >= $a2 && $b3 <= $a3 )
    # # Either the left halfstem overlap or the right halfstem
    # #
    # # USAGE
    # #   $nmerged = mergeCase3( $gap_limit, $nodelist );
    # #   Note: uses global variable $STEM_CENTER_THRESHOLD
    # # ------------------------------------------------------------------------------
    # sub
    # mergeCase4
    # {
    # my( $threshold,
    #
    #
    # @nodelist
    #
    # ) =
    #
    # @_
    #
    # ;
    #
    # my $nmerged = 0;
    # for my $i ( @ nodelist ) {
    # # skip already merged stems and the root node
    # next if ( $i->{id} =~ / x / | | $i->{id} == 0 );  # what does this mean? $i->{id} == 0
    # # the parent stem
    # my ( $a0, $a1, $a2, $a3, $ac ) = ( $i->{pos}->[0], $i->{pos}->[1],
    # $i->{pos}->[2], $i->{pos}->[3],
    # $i->{pos}->[4]
    # );
    #
    # my @ children = @ {$i->{children}};
    # foreach my $child ( @ children ) {
    # next if ( $child->{id} =~ / x / );
    # my ( $b0, $b1, $b2, $b3, $bc ) = ( $child->{pos}->[0], $child->{pos}->[1],
    # $child->{pos}->[2], $child->{pos}->[3],
    # $child->{pos}->[4]
    # );
    #
    # # The left halfstem overlaps, but the right does not, implicit b1< a2
    # if ($b0 <= $a1 & & $b0 >= $a0 & & $b3 < $a2 & & $b2 > $a1 ) {
    # my $stemcenter_diff = abs($bc - $ac);
    # my $right_halfstem_gap = $a2 - $b3;
    # if ( $right_halfstem_gap <= $threshold and
    # $stemcenter_diff <= $STEM_CENTER_THRESHOLD ) {
    #
    # printf STDERR "   Case 4: Child node %s [ %d %d %d %d ] merged into ",
    # $child->{id}, $b0, $b1, $b2, $b3 unless $quiet;
    #
    # printf STDERR "%s [ %d %d %d %d ]\n",
    # $i->{id}, $a0, $a1, $a2, $a3, $ac  unless $quiet;
    #
    # treeMergeNodes( $i, $child );
    # $nmerged++;
    # }
    # }
    #
    # # The right halfstem overlaps, but the left does not, implicit b2> a1
    # if ($b0 > $a1 & & $b1 < $a2 & & $b3 >= $a2 & & $b3 <= $a3) {
    # my $stemcenter_diff = abs($bc - $ac);
    # my $left_halfstem_gap = $b0 - $a1;
    # if ( $left_halfstem_gap <= $threshold and
    # $stemcenter_diff <= $STEM_CENTER_THRESHOLD ) {
    #
    # printf STDERR "   Case 4: Child node %s [ %d %d %d %d ] merged into ",
    # $child->{id}, $b0, $b1, $b2, $b3 unless $quiet;
    #
    # printf STDERR "%s [ %d %d %d %d ]\n",
    # $i->{id}, $a0, $a1, $a2, $a3, $ac  unless $quiet;
    #
    # treeMergeNodes( $i, $child );
    # $nmerged++;
    # }
    # }
    # }
    # }
    # return $nmerged;
    # }
    #
    # # end of mergeCase4
    #
    # # ------------------------------------------------------------------------------
    # # mergeCase5
    # #
    # #       |0aaaaaaaa1|                    |2aaaaaaaa3|
    # #     |0bbbb1|                            |2bbbb3|
    # #
    # #       |0aaaaaaaa1|                    |2aaaaaaaa3|
    # #          |0bbbb1|                           |2bbbb3|
    #
    # #
    # # Looks for stems that overlap and are both nested within the same parent.  For
    # # this case a and b are children of the same parent, d.
    # #
    # # merge non-nested stems, when both halfstem overlaps
    # # number of matched and unmatched bases are below a threshold
    # #
    # # The stem center difference also has to be below a threshold
    # #
    # # ( $b0 < $a0 && $b1 >= $a0 && $b3 <= $a3 && $b3 >= $a2 ) ||
    # # ( $b3 > $a3 && $b2 <= $a3 && $b0 >= $a0 && $b0 <= $a1 )
    # # plus, satisfy, matched, unmatched and stem center threshold
    # #
    # # USAGE
    # #   $nmerged = mergeCase5( $threshold, $nodelist );
    # # ------------------------------------------------------------------------------
    # sub
    # mergeCase5
    # {
    # my( $threshold,
    #
    #
    # @nodelist
    #
    # ) =
    #
    # @_
    #
    # ;
    #
    # my $nmerged = 0;
    # for my $i ( @ nodelist ) {
    # # skip already merged stems
    # next if ( $i->{id} =~ / x / );
    #
    # # parent stem
    # my ( $d0, $d1, $d2, $d3, $dc ) = ( $i->{pos}->[0], $i->{pos}->[1],
    # $i->{pos}->[2], $i->{pos}->[3],
    # $i->{pos}->[4]
    # );
    #
    # my @ children = @ {$i->{children}};
    # foreach my $c_one ( 0..$  # children ) {
    # my $c1 = $children[$c_one];
    # next if ( $c1->{id} =~ / x / );  # skip already merged stems
    #
    # # child 1 stem
    # my ( $a0, $a1, $a2, $a3, $ac ) = ( $c1->{pos}->[0], $c1->{pos}->[1],
    # $c1->{pos}->[2], $c1->{pos}->[3],
    # $c1->{pos}->[4]
    # );
    #
    # my ($left_unmatched, $left_matched)   = (0, 0);
    # my ($right_unmatched, $right_matched) = (0, 0);
    #
    # foreach my $c_two ( $c_one+1..$  # children ) {
    # my $c2 = $children[$c_two];
    # next if ( $c2->{id} =~ / x / );  # skip already merged stems
    #
    # # child 2 stem
    # my ( $b0, $b1, $b2, $b3, $bc ) = ( $c2->{pos}->[0], $c2->{pos}->[1],
    # $c2->{pos}->[2], $c2->{pos}->[3],
    # $c2->{pos}->[4]
    # );
    #
    # #
    # #
    # if ( ( $b0 < $a0 & & $b1 >= $a0 & & $b3 <= $a3 & & $b3 >= $a2 ) | |
    # ( $b3 > $a3 & & $b2 <= $a3 & & $b0 >= $a0 & & $b0 <= $a1 )    ) {
    #
    # # _in and _out are inner and outer boundries of the halfstem
    # # there are two possible cases depending on whether a is to
    # # the left or to the right of b
    #
    # my $l1_in = $b0; my $l1_out = $a0;  # left halfstem, left inner and outer boundries of the halfstem
    # if ( $a0 > $b0 ) {$l1_in = $a0; $l1_out = $b0;}
    #
    # my $l2_in = $b1; my $l2_out = $a1;  # left halfstem, right inner and outer boundries of the halfstem
    # if ( $a1 < $b1 ) {$l2_in = $a1; $l2_out = $b1;}
    #
    # my $r1_in = $b2; my $r1_out = $a2;  # right halfstem, left inner and outer boundries of the halfstem
    # if ( $a2 > $b2 ) {$r1_in = $a2; $r1_out = $b2;}
    #
    # my $r2_in = $b3; my $r2_out = $a3;  # right halfstem, right inner and outer boundries of the halfstem
    # if ( $a3 < $b3 ) {$r2_in = $a3; $r2_out = $b3;}
    #
    # # different length stems would create a imbalance in the
    # # unmatched base count, so it needs to be normalized
    # my $left_stem_len_diff  = abs(($a1-$a0) - ($b1-$b0));
    # my $right_stem_len_diff = abs(($a3-$a2) - ($b3-$b2));
    # $left_unmatched   = $l1_in - $l1_out + $l2_out - $l2_in;
    # $left_unmatched -= $left_stem_len_diff;
    #
    # $right_unmatched  = $r1_in - $r1_out + $r2_out - $r2_in,
    # $right_unmatched -= $right_stem_len_diff;
    #
    # my $c1_len = int(($a1-$a0 + 1 + $a3 - $a2 + 1) / 2);  # average of stem length
    # my $c2_len = int(($b1-$b0 + 1 + $b3 - $b2 + 1) / 2);
    # my $min_stem_len = 0;
    # $min_stem_len = $c2_len;
    # if ( $c1_len < $c2_len ) {$min_stem_len = $c1_len;}  # smaller of the two average of stem length
    #
    # $left_matched  = $l2_in - $l1_in + 1;  # atleast two bases have to be matched
    # $right_matched = $r2_in - $r1_in + 1;
    # # it can also determined from stem length or 30% (which ever
    # # is larger). the unmatched parameter become a bit tricky,
    # # when merging a long stem with a smaller stem.
    # my $MIN_MATCH_THRESHOLD = 1;  # for stem length of 3,2
    # if ( $min_stem_len > 3 ) {
    # $MIN_MATCH_THRESHOLD = 2;  # it remains 2, for 4,5,6
    # if ( 0.3 * $min_stem_len > 2 ) {
    # $MIN_MATCH_THRESHOLD = int( 0.3 * $min_stem_len );  # 1/3 of average stem length
    # }
    # }
    #
    # # Unmatched threshold can be 6
    # # unmatch threshold can be 50% of half stem length or which
    # # ever is smaller... Times 2 or 8
    # $threshold = $min_stem_len;  # 50% of base pairs can be unmatched
    #
    # my $stemcenter_diff = abs($bc - $ac);
    #
    # if ( $left_unmatched <= $threshold and
    # $left_matched >= $MIN_MATCH_THRESHOLD and
    # $right_unmatched <= $threshold and
    # $right_matched >= $MIN_MATCH_THRESHOLD and
    # $stemcenter_diff <= $STEM_CENTER_THRESHOLD + 1 ) {
    #
    # printf STDERR "   Case 5: Child node %s [ %d %d %d %d ]merged to",
    # $c1->{id}, $a0, $a1, $a2, $a3, $ac  unless $quiet;
    # printf STDERR "child node %s [ %d %d %d %d ]\n",
    # $c2->{id}, $b0, $b1, $b2, $b3 unless $quiet;
    #
    # treeMergeNodes( $c1, $c2 );
    # $nmerged++;
    # }
    # }
    # }
    # }
    # }
    # return $nmerged;
    # }
    #
    # # end of mergeCase5

    @staticmethod
    def sample(adj, n):
        """-----------------------------------------------------------------------------------------
        randomly sample a connected subgraph of size=n from the topology.  The sampled graph will be
        connected by non-s (and currently non-x) edges, but size may be less than n if the graph
        being sampled is smaller than size or has disconnected segments. Because there 9 graphs
        with three or fewer stems, it makes little sense to go lower than 3. Sampling is based on
        the adjacency matrix.

        :param n: int   size of graph to sample
        :param n: int   minimum size for sampled graph
        :return: list   list of vertices in sampled graph
        -----------------------------------------------------------------------------------------"""
        attempt_max = 50
        nvertex = len(adj)
        random.seed()

        vlist = []
        excluded = set()
        # initialize neighbor with a random vertex for the first vertex
        # neighbor = set([random.randrange(nvertex)])
        neighbor = {random.randrange(nvertex)}
        attempt = 0
        while len(vlist) < n:
            # print(attempt)
            attempt += 1
            # randomly determine starting vertex
            v0 = random.sample(neighbor, 1)[0]
            vlist.append(v0)
            neighbor.discard(v0)
            excluded.add(v0)
            # update neighbor and exclude sets
            for a in range(nvertex):
                if adj[v0][a] == 'x':
                    neighbor.discard(a)
                    excluded.add(a)
                if adj[v0][a] in 'ijo':
                    neighbor.add(a)

            # make sure no excluded are in the neighbor set
            neighbor = neighbor - excluded

            if len(vlist) >= n:
                # sampled graph is large enough, exit sampling
                break

            if not neighbor:
                # there are no more neighbors but the graph is too small
                # reinitialize vlist, neighbor, and excluded and try again

                vlist = []
                neighbor = {random.randrange(nvertex)}
                excluded.clear()
                if attempt > attempt_max:
                    break

        # end of size < n loop

        return sorted(vlist)

    @staticmethod
    def sample2(adj, n):
        """-----------------------------------------------------------------------------------------
        randomly sample a connected subgraph of size=n from the topology.  The sampled graph will be
        connected by non-s (and currently non-x) edges, but size may be less than n if the graph
        being sampled is smaller than size or has disconnected segments. Because there 9 graphs
        with three or fewer stems, it makes little sense to go lower than 3. Sampling is based on
        the adjacency matrix.

        :param n: int    size of graph to sample
        :param n: int   minimum size for sampled graph
        :return: list   list of vertices in sampled graph
        -----------------------------------------------------------------------------------------"""
        nvertex = len(adj)

        random.seed()
        size = 0

        too_few_nbor_ct = 0
        vlist = []
        neighbor = [random.randrange(nvertex)]
        size = 0
        while size < n:
            # randomly determine starting vertex
            v0 = random.choice(neighbor)
            vlist.append(v0)
            neighbor.remove(v0)
            # print(f'topology:sample v0={v0}\tnvertex={nvertex}')
            size += 1

            for a in neighbor:
                if adj[v0][a] == 'x':
                    neighbor.remove(a)

            # update list of neighbors from adjacency matrix
            v0adj = adj[v0]
            for a in range(len(v0adj)):
                if v0adj[a] in 'sx0-':
                    # skip s and x edges, and self
                    continue

                if (a in neighbor) or (a in vlist):
                    # skip if already in neighbor or vlist
                    continue

                is_x = False
                for v in vlist:
                    # exclude edge if it is exclusive with any previous edges
                    # because otherwise you can end up with graphs that are not ioj connected
                    if adj[v][a] == 'x':
                        is_x = True
                if is_x:
                    continue
                else:
                    # passed all tests, add a to neighbor list
                    neighbor.append(a)

            # if not neighbor:
            #     # neighbor list is empty
            #     break

            if len(vlist) >= n:
                break

            if not neighbor:
                # if there are no more neighbors, you must stop and start again, the outer
                # loop checks to make sure the sample graph is at least size == min_n
                # v0 = random.randrange(nvertex)
                # print(f'v0={v0}\tnvertex={nvertex}')
                # print(f'topology:sample retry\tsize:{size}\tneighbor:{neighbor}\tvlist:{vlist}\tn:{n}')
                vlist = []
                neighbor = [random.randrange(nvertex)]
                size = 0

                too_few_nbor_ct += 1
                if too_few_nbor_ct < 20:
                    continue
                else:
                    break

            # select new vertex and remove from current neighbor list
            # v0 = random.choice(neighbor)
            # neighbor.remove(v0)

            # end of size < n loop

        # end of test for size >= n

        vlist.sort()
        return vlist

    @staticmethod
    def samplebyweight(adj, n, w):
        """-----------------------------------------------------------------------------------------
        randomly sample a connected subgraph of size=n from the topology weighted by w.  The sampled
        graph will be
        connected by non-s edges, but size may be less than n if the graph being sampled is smaller
        than size or has disconnected segmentsSampling is based on the adjacency matrix.

        :param n:
        :return: topology
        -----------------------------------------------------------------------------------------"""
        # randomly determine starting vertex
        # adj = self.adjacency
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
            v0 = random.choices(neighbor, weights=ww)[0]
            neighbor.remove(v0)

        # end of vertex selection loop
        vlist.sort()

        return vlist

    def sample_topology(self, n):
        """-----------------------------------------------------------------------------------------
        Return a new topology sampled from the current one.
        :param n: int, number of stems to sample
        :return: Topology
        -----------------------------------------------------------------------------------------"""
        vlist = self.sample(self.adjacency, n)

        # build the new topology from the current one
        newtopo = Topology()
        for s in self.stem_list:
            if int(getattr(s, 'name')) in vlist:
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

    def sample_xios(self, n, retry=10):
        """-----------------------------------------------------------------------------------------
        Return a xios structure sampled from the current topology. From the selected vertices
        find all the edges that connect those vertices. If the graph is small or disjoint, the
        sampled vlist can be empty.  Try up to retry times and then give up.

        :param self: Xios object
        :param n: int, number of stems to sample
        :param retry: int, number of times to retry if vlist is too small
        :return: Xios object (see xios.py)
        -----------------------------------------------------------------------------------------"""
        from xios import Xios

        edge = {'i': 0, 'j': 1, 'o': 2, 's': 3, 'x': 4}

        vlist = []
        tries = 0
        while len(vlist) < n and tries < retry:
            tries += 1
            vlist = Topology.sample(self.adjacency, n)

            adj = self.adjacency
            struct = []

            # identify all the edges between the vertices in vlist
            for r in range(len(vlist) - 1):
                row = vlist[r]
                for c in range(r + 1, len(vlist)):
                    col = vlist[c]
                    if adj[row][col] in 'ijo':
                        struct.append([row, col, edge[adj[row][col]]])

        # print(struct)
        return Xios(list=struct)

    def sample_xios_weighted(self, n, w):
        """-----------------------------------------------------------------------------------------
        Return a xios structure sampled from the current topology.

        :param n: int, number of stems to sample
        :param w: list of weights
        :return: Xios object (see xios.py)
        -----------------------------------------------------------------------------------------"""
        from xios import Xios

        adj = self.adjacency
        edge = {'i': 0, 'j': 1, 'o': 2, 's': 3, 'x': 4}
        vlist = self.samplebyweight(adj, n, w)


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

    def __init__(self, list=None):
        """-----------------------------------------------------------------------------------------
        SerialRNA is a list with some extra methods.  Technically it doesn't need a constructor
        since it inherits one from list, but just in case attributes need to be added, this
        constructor just calls the parent's
        -----------------------------------------------------------------------------------------"""
        if list is None:
            list = []
        super(SerialRNA, self).__init__(list)


        # add additional attributes

    def connected(self):
        """-----------------------------------------------------------------------------------------
        return the connected graph(s) in the current structure. The entire graph is connected if one
        structure is returned

        :return: list of SerialRNA
        -----------------------------------------------------------------------------------------"""
        component = []
        openstem = []
        begin = 0
        for pos in range(len(self)):
            if self[pos] in openstem:
                openstem.remove(self[pos])
                if len(openstem) == 0:
                    component.append(SerialRNA(self[begin:pos + 1]))
                    begin = pos + 1
            else:
                openstem.append(self[pos])

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
                # begin and end are the locations of the new stem (stem 0)
                extended_rna = [None for _ in range(len(self) + 2)]
                extended_rna[begin] = 0
                extended_rna[end] = 0
                newpos = 0
                for pos in self:
                    # copy in the old stems, except where the new stem is placed
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
        convert = {}
        for pos in self:
            if pos not in convert:
                convert[pos] = stem
                stem += 1

        for pos in range(len(self)):
            self[pos] = convert[self[pos]]

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
    think of a topology with n stems linear set of coordinates from 0 to n-1 pair format lists the
    coordinates or each stem in this abstract coordinate system.  For instance, the abstract
    structure (()()) would be (0,5) (1,2) (3,4) or the pseudoknotted structure {[(])) would be
    (0,5) (1,3) (2,4).  This representation makes it very simple to add new stems to an existing
    structure.

    Internally stored as a simple list instead of a list of lists (pairs).  so the above stuctures
    are simply [0,5,1,2,3,4] and [0,5,1,3,2,4]

    Synopsis
        from graph import PairRNA

        graph = PairRNA()
        graph.fromList([0,1,1,0])
            or
        graph = PairRNA(l=[0,1,1,0])

        graph.reverse()     # reverse the order of stems left-right
    ============================================================================================="""

    def __init__(self, topology=None, topology_type='pair'):
        """-----------------------------------------------------------------------------------------
        Internal data structure is a list of lists:
        each element in the list of the begin, end position of a stem.  The stems
        are the coordinates of the left and right half stems as integers along a line.

        two nested stems are [ [0,3], [1,2] ]
        a simple pseudoknot is [ [0,2], [1,3] ]

        :param inlist: a list of lists with coordinates of left and right half stems
        -----------------------------------------------------------------------------------------"""
        if topology is None:
            # avoid mutable initialization in arguments
            topology = []
        self.pairs = []
        self.nstem = 0

        if topology:
            if topology_type == 'pair':
                self.pairs = topology
            elif topology_type == 'serial':
                self.from_SerialRNA(topology)
            else:
                sys.stderr.write(f'PairRNA.init - unknown topology type {topology_type}\n')

    def __str__(self, sep='_'):
        """-----------------------------------------------------------------------------------------
        Return a serialized version of the pair structure
        :return: string
        -----------------------------------------------------------------------------------------"""
        s = ''
        for i in range(0, len(self.pairs), 2):
            try:
                s += '({},{}) '.format(self.pairs[i], self.pairs[i + 1])
            except IndexError:
                # undefined stem
                s += '(.,.) '

        return s.rstrip('_')

    def __len__(self):
        """-----------------------------------------------------------------------------------------
        length is the number of stems
        -----------------------------------------------------------------------------------------"""
        return len(self.pairs) // 2

    def from_SerialRNA(self, serial):
        """"----------------------------------------------------------------------------------------
        read in a graph as a SerialRNA object.  A Serial RNA list will look like [0,1,1,2,2,0]
        which is the structure (()()) and the PairRNA [0,5,1,2,3,4]

        :param serial: SerialRNA object, structure in list format
        :return: int, number of stems
        -----------------------------------------------------------------------------------------"""
        self.nstem = int(len(serial) // 2)
        self.pairs = [0 for _ in range(len(serial))]

        # make sure the graph is canonically numbered
        next = 0
        seen = {}
        for i in range(len(serial)):
            pos = serial[i]
            if pos not in seen:
                seen[pos] = next
                next += 1
            serial[i] = seen[pos]

        # the values indicate the stem number, the index indicates the position
        pos = [p * 2 for p in range(self.nstem)]
        for i in range(len(serial)):
            stem = serial[i]
            self.pairs[pos[stem]] = i
            pos[stem] += 1

        return self.nstem

    def from_SerialRNA_string(self, serialstr, sep=' '):
        """-----------------------------------------------------------------------------------------
        the input list is a string separated by sep, e.g. '0 1 1 0'

        :param serialstr: string, input graph as a SerialRNA string
        :param sep: string, separation character in input string
        :return: integer, number of stems
        -----------------------------------------------------------------------------------------"""
        pairs = self.pairs
        serial = []
        for s in serialstr.split(sep):
            serial.append(int(s))

        # make sure the graph is canonically numbered
        next = 0
        seen = {}
        for i in range(len(serial)):
            pos = serial[i]
            if pos not in seen:
                seen[pos] = next
                next += 1
            serial[i] = seen[pos]

        self.pairs = [0 for p in range(len(serial))]
        pos = [p * 2 for p in range(self.nstem)]
        for i in range(len(serial)):
            stem = int(serial[i])
            self.pairs[pos[stem]] = i
            pos[stem] += 1

        return len(pairs)

    def to_SerialRNA(self):
        """-----------------------------------------------------------------------------------------
        return the structure in SerialRNA format

        :return g: SerialRNA, structure in list format
        -----------------------------------------------------------------------------------------"""
        g = [0 for _ in range(self.nstem * 2)]

        stem = 0
        for pair in self.pairs:
            g[pair[0]] = stem
            g[pair[1]] = stem
            stem += 1

        return SerialRNA(g)

    def from_vienna(self, vienna):
        """-----------------------------------------------------------------------------------------
        Read structure in Vienna format, e.g., (()) or ([)]. each stem should be represented as a
        single pair of brackets so this is basically the abstract shapes format with pseudoknots.
        Dots are optional, they will produce a struture with non-consecutive positions which can be
        fixed with canonical()

        :param vienna: string, structure in Vienna dot-bracket format
        :return: integer, number of stems
        -----------------------------------------------------------------------------------------"""
        left = {'(': ')', '[': ']', '{': '}', '<': '>'}
        right = ')]}>'

        # find length by count opening brackets
        nstem = 0
        for bracket in vienna:
            if bracket in left:
                nstem += 1
        pairs = [0] * nstem * 2

        stack = {ch: [] for ch in right}
        level = []
        stem = 0
        pos = 0
        for bracket in vienna:
            if bracket in left:
                # opening bracket
                pairs[stem * 2] = pos
                stack[left[bracket]].append(stem)
                stem += 1
            else:
                s = stack[bracket].pop()
                pairs[s * 2 + 1] = pos

            pos += 1

        self.pairs = pairs
        self.nstem = len(self.pairs) // 2

        return self.nstem

    def to_vienna(self, pad=''):
        """-----------------------------------------------------------------------------------------
        return a string with the structure in vienna format.  This is basically the abstract shapes
        format with support for pseudoknots.

        :param pad: string, separator to put between brackets
        :return: string
        -----------------------------------------------------------------------------------------"""
        level = {'(': 0, ')': 0, '[': 1, ']': 1, '{': 2, '}': 2, '<': 3, '>': 3}
        left = '([{<'
        right = ')]}>'
        nbracket = len(left)

        pairs = self.pairs
        vienna = [None for _ in range(self.nstem * 2)]

        for pos in range(0, len(pairs), 2):
            lpos = pairs[pos]
            rpos = pairs[pos + 1]

            count = [0] * nbracket
            for pos in range(lpos, rpos):
                # look at the range of this stem and count the number of each
                # type brackets that are in use
                if vienna[pos]:
                    count[level[vienna[pos]]] += 1

            lbracket = rbracket = ''
            for i in range(nbracket):
                if not count[i] % 2:
                    # if the count is even, this bracket type can be used. otherwise it's a pknot
                    lbracket = left[i]
                    rbracket = right[i]
                    break

            vienna[lpos] = lbracket
            vienna[rpos] = rbracket

        return pad.join(vienna)

    def reverse(self):
        """-----------------------------------------------------------------------------------------
        reverse the positions of the stems by converting to maxpos-pos and resorting in order
        of first coordinate

        :return: PairRNA object, a new object in reversed order
        -----------------------------------------------------------------------------------------"""
        m = len(self.pairs) - 1
        new = self.duplicate()

        for i in range(0, m, 2):
            # rev[0] = m - self.pairs[i][0]
            # rev[1] = m - self.pairs[i][1]
            # self.pairs[i][0], self.pairs[i][1] = rev[1], rev[0]
            new.pairs[i], new.pairs[i + 1] = m - self.pairs[i + 1], m - self.pairs[i]

        new.canonical()

        return new

    def duplicate(self):
        """-----------------------------------------------------------------------------------------
        Make a deep copy of the current topology (self)
        :return: Topology.PairRNA object
        -----------------------------------------------------------------------------------------"""
        new = PairRNA()
        new.pairs = copy.deepcopy(self.pairs)
        new.nstem = len(new.pairs) // 2

        return new

    def canonical(self):
        """-----------------------------------------------------------------------------------------
        reorder the stem pairs based on the first coordinate of the pair

        :return: True
        -----------------------------------------------------------------------------------------"""
        l = len(self.pairs)

        # make numbering sequential
        map = []
        for s in sorted(self.pairs, key=lambda p: p):
            if s not in map:
                map.append(s)

        for i in range(l):
            # renumber sequentially
            self.pairs[i] = map.index(self.pairs[i])
            if i % 2:
                # check that beginning and end are in increasing order
                if self.pairs[i] < self.pairs[i - 1]:
                    self.pairs[i], self.pairs[i - 1] = self.pairs[i - 1], self.pairs[i]

        # sort by first position in pair
        new = map
        n = 0
        for i in sorted(range(0, l, 2), key=lambda p: self.pairs[p]):
            new[n], new[n + 1] = self.pairs[i], self.pairs[i + 1]
            n += 2

        self.pairs = new

        return True

    def push_pair(self, pair):
        """-----------------------------------------------------------------------------------------
        add a pair (two element list) to the beginning of the pairs list

        :param pair: list, two element list describing the paired positions of a stem
        :return: int, number of stems in list
        -----------------------------------------------------------------------------------------"""
        self.pairs = pair + self.pairs
        self.nstem = len(self.pairs) // 2

        return self.nstem

    def connected(self):
        """-----------------------------------------------------------------------------------------
        Returns True if graph is i-o connected.  Assumes that the PairRNA is in canonical form, i.e.,
        that the starting coordinate of the stem always increases for each stem.

        :return: True / False
        -----------------------------------------------------------------------------------------"""
        pairs = self.pairs
        begin = 0
        end = 0
        for i in range(0, len(pairs), 2):
            if pairs[i] <= end:
                end = max(end, pairs[i + 1])
            else:
                # disconnected
                # print(f'\tnot connected {self}')
                return False

        return True

    def depth(self):
        """-----------------------------------------------------------------------------------------
        Return  list with the nesting depth at each position of the graph in list format.  This is
        useful for mountain plots and determining connectivity.  If the level reaches zero before
        the last position, the graph is disconnected.

        :return: list, number of stems each position of the pseudo sequence is contained in
        -----------------------------------------------------------------------------------------"""
        pairs = self.pairs
        l = len(self.pairs)
        depth = [0 for _ in range(l)]

        for pos in range(0, l, 2):
            for p in range(pairs[pos], pairs[pos + 1] + 1):
                depth[p] += 1

        return depth


####################################################################################################
####################################################################################################
class RNAstructure(Topology):
    """=============================================================================================
    RNA structure class is for working with output from the RNAStructure package
    a single RNA structure
    ============================================================================================="""

    def __init__(self, *args, **kwargs):

        super(RNAstructure, self).__init__(*args, **kwargs)
        self.energy = 0.0  # not defined for probknot
        self.pair = []  # base number of the paired base

    def is_ctheader(self, field):
        """-----------------------------------------------------------------------------------------
        Read and parse the header line for a single structure in a CT file.
        Format differs with the program creating the file.

          240  ENERGY = -49.3   mr_s129                                     <= mfold
          371     dG = -80.1      Mycpl_arthr 371 bases                     <= unafold/mfold
          229  ENERGY = -102.8  rnasep_m.A_fulgidus RNase P RNA length=229  <= RNAstructure/Fold
          240   mrs129                                                      <= RNAstructure/probknot

          Assume the first non-blank line is a header line.

        :param field: list, tokens split from current line of file
        :return: boolean, True if header is parsed
        -----------------------------------------------------------------------------------------"""
        if len(field) == 2:
            # probknot
            self.sequence_length = int(field[0])
            self.sequence_id = field[1]
            return True

        if field[1] == 'ENERGY' or field[1] == 'dG':
            try:
                self.energy = float(field[3])
                self.sequence_id = field[4]
                self.sequence_length = int(field[0])
                return True
            except:
                line = '\t'.join(field)
                sys.stderr.write(
                    f'Topology/RNAstructure::is_ctheader - cannot read header ({line})\n')

        return False

    def is_ctdata(self, field):
        """-------------------------------------------------------------------------------------
        A CT data row has exactly six fields: fields 0 and 2-5 are integers, field 1 is a
        single character

        :param field: list of strings - fields in original line
        :return: boolean
        -------------------------------------------------------------------------------------"""
        ok = True
        if not len(field) == 6:
            ok = False

        if not field[1].isalpha() and len(field[1]) == 1:
            ok = False

        else:
            for i in (0, 2, 3, 4, 5):
                try:
                    field[i] = int(field[i])
                except ValueError:
                    ok = False

        if not ok:
            content = '\t'.join(field)
            sys.stderr.write(
                f'Topology/RNAstructure::is_ctdata - error in dataline ({content})\n')
            exit(1)

        return ok

    def CTRead(self, filename, ddG=4.0):
        """-----------------------------------------------------------------------------------------
        Read in RNAstructure CT file

        format example (unifold/mfold):
          240  ENERGY = -49.3   mr_s129
            1 A       0    2    0    1
            2 A       1    3    0    2
            3 C       2    4    0    3
            4 C       3    5    0    4
            5 A       4    6    0    5

            base_number base previous_base next_base paired_base base_number
        format example (probknot)
          240   mrs129
            1 A       0    2    0    1
            2 A       1    3    0    2
            3 C       2    4    0    3
            4 C       3    5    0    4

        CT files written by Fold in the RNAstructure package have the heading
           229  ENERGY = -102.8  rnasep_m.A_fulgidus RNase P RNA length=229

        usage
            rna.CTRead(filename)

        :param filename: string, filename with CT formatted structure
        :param ddG: float, maximum delta delta-G from mfe to include
        :return: integer, number of bases
        -----------------------------------------------------------------------------------------"""
        ct = None
        try:
            ct = open(filename, 'r')
        except OSError:
            sys.stderr.write("Topology/RNAstructure::CTRead - unable to open CT file ({})".format(filename))

        nbase = 0
        path = filename.split('/')
        self.filename = path[-1]

        nstruct = 0
        comment = f'\tTopology/RNAStructure::CTRead - Converting CT file {filename} to XIOS delta(deltaG)={ddG}\n'
        mfe = 0.0
        for line in ct:
            if not line:
                # skip blank lines?
                continue
            field = line.split()

            if self.is_ctheader(field):
                # add provenance information to Topology comment
                if nstruct == 0:
                    # for the first structure the pair array is empty, save mfe and create pair list
                    # pairlist is the workspace for the stem calculation
                    mfe = float(self.energy)
                    self.pair = [0] * (self.sequence_length + 1)

                else:
                    # second or later structure, save the stems, check that we are below mfe + ddG
                    # 0.01 is a fudge factor to make sure floating point differences don't
                    # cause problems
                    if self.energy > mfe + ddG + 0.01:
                        break
                    self.stemlist_from_pairs(unpaired=2)
                    self.pair = [0] * (self.sequence_length + 1)

                comment += f'\tstructure {nstruct}: {line}\n'
                nstruct += 1


            elif self.is_ctdata(field):
                base = field[1]
                n = int(field[0])
                pair = int(field[4])
                self.sequence += base
                if pair != 0:
                    self.pair[pair] = n
                    self.pair[n] = pair
                nbase += 1

            else:
                # unidentified line
                sys.stderr.write(f'Topology/RNAstructure::CTread - unidentified line ({line})\n')

        # add the final structure, energy is checked vs ddG before reading pairs, above
        self.stemlist_from_pairs(unpaired=2)
        dedup_n = self.stemlist_deduplicate()
        comment += f'\t{dedup_n} stems after removing duplicates'
        self.comment.append(comment)

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

    def stemlist_deduplicate(self):
        """-----------------------------------------------------------------------------------------
        Remove stems from the stemlist that are identical or perfectly nested.  The nesting test is
        based only on ranges so the base-pairing may differ slightly

        TODO: move to Topology since this is a generic function useful for any Topology

        :return: int, number of stems in revised stemlist
        -----------------------------------------------------------------------------------------"""
        s = self.stem_list
        new = []

        order = sorted(range(len(s)), key=lambda x: (len(s[x].lvienna)), reverse=True)
        order = sorted(order, key=lambda x: (s[x].lbegin, s[x].rbegin))
        # for i in order:
        #     print(s[i].lbegin, s[i].lend, s[i].rbegin, s[i].rend)

        skip = [False] * len(s)
        i = 0
        j = 1
        while i < len(order) - 1:
            if skip[i]:
                i += 1
                j = i + 1
                continue

            s1 = s[order[i]]
            s2 = s[order[j]]
            # if s2.lbegin > s1.lend:
            #     i += 1
            #     continue

            while s2.lbegin <= s1.lend and j < len(order):
                # check potential overlap
                if s2.lend <= s1.lend and s2.rbegin >= s1.rbegin and s2.rend <= s1.rend:
                    skip[j] = True
                    j += 1
                    if j < len(order):
                        s2 = s[order[j]]
                else:
                    break

            i += 1
            j = i + 1

        nstem = 0
        for i in range(len(order)):
            stem = s[order[i]]
            if not skip[i]:
                stem.name = nstem
                new.append(stem)
                # print(stem.lbegin, stem.lend, stem.rbegin, stem.rend, stem.lvienna)
                nstem += 1

        self.stem_list = new
        return nstem

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
                # is this position unpaired (0), or the right half stem
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
            stem.name = '{}'.format(nstem)
            self.stem_list.append(stem)
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
            self.stem_list.append(s)
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
        for stem in self.stem_list:
            n += 1
            stemstr += '{0}\t{1}\n'.format(n, stem.formatted())
        return stemstr

    def adjacency_from_stemlist(self):
        """-----------------------------------------------------------------------------------------
        Calculate an adjacency matrix from stem_list. assumes stems are ordered by the beginning
        of the left half-stem.
        :return: dict, keys are edge types, values are counts
        -----------------------------------------------------------------------------------------"""
        nstem = len(self.stem_list)

        edges = {'i': 0, 'j': 0, 'o': 0, 's': 0, 'x': 0}
        a = [[0 for _ in range(nstem)] for _ in range(nstem)]

        for i in range(nstem):
            stem_i = self.stem_list[i]

            for j in range(i + 1, nstem):
                stem_j = self.stem_list[j]

                if stem_i.rend < stem_j.lbegin:
                    # serial edge
                    a[i][j] = 's'
                    a[j][i] = 's'
                    edges['s'] += 1

                # elif stem_i.lend < stem_j.lbegin and stem_i.rend < stem_j.rbegin:
                elif stem_i.lend < stem_j.lbegin and stem_i.rbegin > stem_j.lend and stem_i.rend < stem_j.rbegin:
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
        for edge in self.edge_list:
            # include and whole arguments don't do anything
            n += 1
            edgestr += f'{n}: '
            for neighbor in edge:
                edgestr += '{neighbor[0]}{neighbor[1]}, '

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
        pairgraphs = [[2, 3, 0, 4], [0, 3, 1, 2], [0, 2, 1, 4, 3, 5], [0, 4, 2, 3], [1, 0, 3, 5]]
        serialgraphs = [[0, 0, 1, 1, 2, 2], [0, 1, 0, 1, 2, 2], [3, 1, 1, 3], [3, 2, 0, 2, 0, 3]]
        viennagraphs = ['(())', '()(())', '([)]', '<([>])', '<[(>])', '([)(])']

        # PairRNA
        print('Testing PairRNA')
        print('\nPairRNA graphs')
        for g in pairgraphs:
            print(f'Pair {g} =>', end=' ')
            pair = PairRNA()
            pair.pairs = g
            pair.canonical()
            print(f'{pair.pairs}')

        print('\nSerialRNA graphs')
        for g in serialgraphs:
            pair.from_SerialRNA(g)
            print('Serial {} => {}'.format(g, pair))

        gstring = '0,2,2,0'
        pair.from_SerialRNA_string(gstring, sep=',')
        print('Serial string {} => {}'.format(gstring, pair))

        print('\nVienna format')
        for v in viennagraphs:
            pair.from_vienna(v)
            print(f'Vienna {v} => {pair.pairs}', end='')
            vienna = pair.to_vienna()
            print(f' => {vienna}   depth={pair.depth()}')

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


    def test_SerialRNA():
        """-----------------------------------------------------------------------------------------
        SerialRNA format is a list showing the begin and end points of each numbered stem in an
        abstract sequentially numbered sequence. Examples
        (()) = [0,1,1,0]
        ([)] = [0,1,0,1]
        ([)(]) = [0,1,0,2,1,2]
        -----------------------------------------------------------------------------------------"""
        rnas = [[0, 0, 1, 1, 2, 2],
                [0, 1, 0, 1, 2, 2],
                [0, 1, 1, 2, 2, 0],
                [0, 1, 2, 1, 2, 0],
                [0, 0], []
                ]

        print('canonical form')
        noncanonical = [[3, 3, 0, 0, 1, 1], [1, 1, 2, 2, 3, 3], [1, 1, 2, 2, 4, 4],
                        [3, 2, 0, 2, 0, 3]]
        empty = SerialRNA()

        for testcase in noncanonical:
            rna = SerialRNA(testcase)
            print('RNA {}'.format(rna))
            rna.canonical()
            print('\tcanonical {}'.format(rna))

        print('\nConnected components')
        for testcase in rnas:
            rna = SerialRNA(testcase)
            print('RNA {}'.format(rna))
            connected = rna.connected()
            if len(connected) > 1:
                for i in range(len(connected)):
                    print('\tcomponent {}: {}'.format(i, connected[i]))

        print('\nExtension')
        for testcase in rnas:
            rna = SerialRNA(testcase)
            print('RNA {}'.format(rna))
            for new in rna.addstemzero():
                print('\t{} {}'.format(new, len(new.connected())))

        return


    def test_RNAstructure(filename):
        """-----------------------------------------------------------------------------------------

        :param filename:
        :return:
        -----------------------------------------------------------------------------------------"""
        test = RNAstructure()
        ddG = 1.75
        test.CTRead(filename, ddG)
        test.XIOSwrite(sys.stdout)

        return


    # ##############################################################################################
    # Testing
    # ##############################################################################################
    test_pair()
    test_topology()
    test_SerialRNA()
    test_RNAstructure('data/mr_s129.fold.ct')

    exit(0)

    # TODO convert this section, or delete
    # for rna in graphs:
    #     g = RNAGraph(rna)
    #     three = RNAstructure()
    #     three.stem_listGetFromGraph(g)
    #     print('stem_list\n')
    #     print(three.stem_list_format())
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
    #     rna.stem_list_from_pairs()
    #     # print(rna)
    #     print('stem_list\n')
    #     print(rna.stem_list_format())
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
