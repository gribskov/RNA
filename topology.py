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


# --------------------------------------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    top = Topology()
    top.XIOSread('data/rnasep_a1.Buchnera_APS.xios')

    exit(0)

"""

#-----------------------------------------------------------------------------
# initialize
#
# This is called by the constructor "new" to set up empty attributes of the
# class.
#
# USAGE
#      $Topology_obj->initialize();
#
#-----------------------------------------------------------------------------

sub initialize {
    my ($self) = @_;
    $self->{stem_list}    = [];
    $self->{vienna}       = '';
    $self->{level_array}  = [];
    $self->{brackets}     = ();
    $self->{sequence}     = '';
    $self->{sequence_id}  = '';
    $self->{sequence_doc} = '';
    $self->{comment}      = [];    # a generic place to store other information
    $self->{length}       = 0;
    return;
}


    foreach my $attr ( keys %$attribute_hash ) {
        if ( $attr =~ /unpaired/i ) {
            $max_unpaired = $$attribute_hash{$attr};

        }
        elsif ( $attr =~ /data|vienna/i ) {
            $self->vienna( $$attribute_hash{$attr} );
            my $length = $self->findStems($max_unpaired);

        }
        elsif ( $attr =~ /mfold|ct/i ) {
            $file = $$attribute_hash{$attr};
            my $n_stem = $self->mfold2( $file, $max_unpaired );

        }
        elsif ( $attr =~ /rnaml/i ) {
            $file = $$attribute_hash{$attr};
            my $n_stem = $self->rnaml( $file, $max_unpaired );

        }
        elsif ( $attr =~ /plot/i ) {
            $file = $$attribute_hash{$attr};
            my $n_stem = $self->energyplot($file);

        }
        elsif ( $attr =~ /pair/ ) {
            my $pairs = $$attribute_hash{pair};
            my $nstem = $self->pairs( $pairs, $max_unpaired );

        }
        elsif ( $attr =~ /parray/ ) {
            my $parray = $$attribute_hash{parray};
            my $nstem  = $self->pairArray($parray);

        }
        else {
            $self->{$attr} = $attribute_hash->{$attr};
        }
    }

    return;
}

        }
    }

    return $found;
}

# end of setBracketTypes

#-----------------------------------------------------------------------------
# addStem
#
# This accessor method adds stem object to the stem list array.
# Output: return;
#
# USAGE
#     $Topology_obj->addStem($stem_obj);
#-----------------------------------------------------------------------------

sub addStem {
    my ( $self, $stem_obj ) = @_;

    push( @{ $self->stemList }, $stem_obj );
    return;
}

#-----------------------------------------------------------------------------
# addTopology
#
# This accessor method adds all the stems in another topology object to the stem list
# array of the existing topology.
#
# USAGE
#     $Topology_obj->addTopology($another_topology_obj);
#-----------------------------------------------------------------------------
sub addTopology {
    my ( $self, $another_topo ) = @_;
    foreach my $stem_obj ( @{ $another_topo->stemList } )  {
        $self->addStem( $stem_obj );
    }
    return;
}
#------------------------------------------------------------------------------
# stemsToPairArray
#
#   Given a list of stems, transfer into pair array format.
#
# USAGE
#   my $pairs = stemsToPairArray( \@stem_list );
#------------------------------------------------------------------------------
sub stemsToPairArray {
    my ($stems) = @_;
    my %coord;
    my @coord;
    my %pairs;
    my @pair = ();
    foreach my $stem ( @{$stems} ) {
        my ( $a1, $a2, $a3, $a4 ) =
          ( $stem->[0], $stem->[1], $stem->[2], $stem->[3] );
        push @coord, ([$a1, $a2], [$a3, $a4] );
        my $ll   = join ".", ( $a1, $a2 );
        my $rr   = join ".", ( $a3, $a4 );
        my $llrr = join ".", ( $ll, $rr );
        $coord{$ll}   = undef;
        $coord{$rr}   = undef;
        $pairs{$llrr} = undef;
    }

    @coord = sort { $a->[0] <=> $b->[0] } @coord;
    #print STDERR "coord: @coord\n";
    my $pos = 0;
    while ( $pos < @coord ) {
        my $cc = join ".", @{$coord[$pos]};
        $coord{$cc} = $pos;
        $pos++;
    }

    foreach my $llrr ( keys %pairs ) {

        my @cc = split /\./, $llrr;
        my $ll = join ".", @cc[ 0 ... 1 ];
        my $rr = join ".", @cc[ 2 ... 3 ];
        my $lc = $coord{$ll};
        my $rc = $coord{$rr};
        #print STDERR "ll $ll rr $rr cc @cc lc $lc rc $rc\n";
        push @pair, ( $lc, $rc );
    }

    return Topology::sortPairs(@pair);

}

# End of stemsToPairArray


    #------------------------------------------------------------------------------
# XIOSWrite
#
# write a XIOS file describing the topology.  This functions is the
# authoritative definition of the format of a XIOS file.  Make sure that any
# changes made here are reflected in readXIOS.
#
# A XIOS file has up to three sections.  Each section is delimited by an XML
# like tag so that it is easier to read them independently.
#
# USAGE
#   $xios_string = $topology->XIOSWrite( );
#   $xios_string = $topology->XIOSWrite( $filename );
#------------------------------------------------------------------------------
sub XIOSWrite {
    my ( $self, $filename ) = @_;
    my $str = "";

    $self->adjacencyMatrix;

    $str .= qq{<XIOS version='$VERSION'>\n};

    $str .= $self->informationFormat;
    $str .= $self->stemlistFormat;
    $str .= $self->edgelistFormat;
    $str .= $self->adjacencyFormat;

    $str .= qq{</XIOS>\n};

    if ( defined $filename ) {

        # if no filname is provided just return the string
        open( my $out, '>', $filename )
          or die "Topology::XIOSWrite unable to open $filename for output\n\n";
        print $out $str;
    }

    return $str;
}

# End of XIOSWrite
#-----------------------------------------------------------------------------
# topologyToVienna
#
# generates dot-bracket (Vienna) representation of RNA structure from
# topology representation (list of stem pairs).
# Takes in number of dots and number of brackets desired in the Vienna
# representation, and the topology array. Returns a string containing
# Vienna representation.
#
# this is used to create a dummy vienna string for a given set of half-stem
# pairs
#
# USAGE:
#   $vienna_rep = topologyToVienna( $ndots, $nbrackets, @topology );
#
# 10 July, 2007
#-----------------------------------------------------------------------------
sub topologyToVienna {
    my ( $self, $ndot, $nbrackets, @topol ) = @_;

    my @out;
    my %bracket = (
        '(' => 0,
        ')' => 0,
        ']' => 1,
        '[' => 1,
        '{' => 2,
        '}' => 2,
        '>' => 3,
        '<' => 3,
    );

    my @pairs = ( [ "(", ")" ], [ "[", "]" ], [ "{", "}" ], [ "<", ">" ], );

    while (@topol) {

        # get the begin and end of a stem
        my $begin = shift @topol;
        my $end   = shift @topol;
        my @used;
        foreach my $i ( ( $begin - 1 ) .. ( $end - 1 ) ) {
            my $level;
            if ( defined $out[$i] ) {

                # if a bracket is present between begin and end of the stem
                # get the level of that bracket from %bracket
                $level = $bracket{ $out[$i] };
                $used[$level]++;
            }
        }
                
        my $found = -1;
        foreach my $level ( 0 .. $#pairs ) {

            # for all 4 levels of brackets in @pairs get the level in use
            if ( !defined( $used[$level] ) || !$used[$level] % 2 ) {
                $found = $level;
                last;
            }
        }
        if ( $found == -1 ) {
            print STDERR "begin=$begin   end=$end   found=$found\n";
            print STDERR "used=@used\n";
            print STDERR "out=@out\n";
            print STDERR "too many bracket levels\n";
            die;
        }
        else {

            # fill in the brackets according to the level in use
            $out[ $begin - 1 ] = $pairs[$found][0];
            $out[ $end - 1 ]   = $pairs[$found][1];
        }
    }

    # join individual brackets to get vienna representation string
    # expand structure by inserting dots and multiples of each bracket
    my $vienna = join( "", @out );
    foreach my $onebrace ( keys %bracket ) {
        my ( $from, $to ) = expandBrackets( $onebrace, $ndot, $nbrackets );
        $vienna =~ s/$from/$to/g;
    }
    return $vienna;
}

      #-----------------------------------------------------------------------------
# mfold
#
# open a file and read the structure from a mfold/ct formatted file.  The ct
# file is converted to a vienna string and then parsed as a vienna string.
# this will not work if the structure contains pseduoknots.
#
# SAMPLE FILE:
#
#369     dG = -138.2     sacce S.cerevisiae RNase P RNA
#1       G       0       2       0       1       0       0
#2       U       1       3       367     2       0       3
#3       G       2       4       366     3       2       4
#4       G       3       5       365     4       3       5
#5       A       4       6       364     5       4       6
#6       A       5       7       363     6       5       7
#
# USAGE
#   $ n_stems  = $topology->mfold( $filename );
#-----------------------------------------------------------------------------
sub mfold2 {
    my ( $self, $filename, $max_unpaired ) = @_;

    open( CT, "<$filename" ) || die "unable to open ct file $filename\n\n";

    my $line = <CT>;
    my ( $length, $dg, $equal, $energy, $id, $doc ) = split " ", $line, 6;
    $self->{energy} = $energy;

    #my @vienna_array;
    #foreach my $i ( 0 .. $length-1 ) {
    #    $vienna_array[ $i ] = ".";
    #}

    # get the basepairs from the CT file

    my @pair;

    my $seq   = "";
    my $npair = 0;
    while ( $line = <CT> ) {

        next if ( $line =~ /#/ );
        last if ( $line =~ /dG/ );    #stop after first structure

                        my ( $left, $base, $pos, $x1, $right, $x2 ) = split " ", $line;
        next if ( ( !( defined($base) ) ) && ( !( defined($pos) ) ) );
        $seq .= $base;

        # in CT file, pairs are listed from both directions
        next unless ( $right && $left < $right );
        $pair[$npair] = [ $left, $right ];
        $npair++;
    }

    #foreach my $ii ( 0..$#pair ) {
    #    print " pair[$ii] = $pair[$ii]->[0]:$pair[$ii]->[1]\n";
    #}
    $self->sequence($seq);
    $self->sequenceId($filename);
    my @stemlist = basepairToStem(@pair);
    @stemlist = splitStems( $max_unpaired, @stemlist );

    my $n_stem = 0;
    foreach my $stem (@stemlist) {
        $self->addStem($stem);

        #print "stem: $n_stem\n";
        #$stem->dump;
        $n_stem++;
    }

    return ($n_stem);
}

#------------------------------------------------------------------------------
# mfoldMultiple
#
# Read the all the suboptimal CT sructures in the specified file down to the
# delta deltaG limit.  If no limit is specified, all structures are read.  The
# list of stems is filtered to remove duplicates.
#
# SAMPLE FILE:
#
#369     dG = -138.2     sacce S.cerevisiae RNase P RNA
#1       G       0       2       0       1       0       0
#2       U       1       3       367     2       0       3
#3       G       2       4       366     3       2       4
#4       G       3       5       365     4       3       5
#5       A       4       6       364     5       4       6
#6       A       5       7       363     6       5       7
#
# USAGE
#    my $n_stem = $topology->mfoldMultiple( $filename, $max_unpaired, $ddG );
#------------------------------------------------------------------------------
sub mfoldMultiple {
    my ( $self, $filename, $max_unpaired, $ddG ) = @_;
    my $n_stem;

    open( my $in, "<", $filename )
      or die
      "Topology::mfoldMultiple, unable to open input CT file ($filename)\n\n";

    $self->sequenceId($filename);

    # read header for first structure

    my $line = <$in>;
    my ( $length, $dg, $equal, $mfe, $id, $doc ) = split " ", $line, 6;
    $self->length($length);

    my @stemlist;
    my ( $seq, @pair );
    while ( my $line = <$in> ) {
        next if ( $line =~ /#/ );    # skip comments
        next unless ( $line =~ /\w/ );    # skip blank lines

        if ( $line =~ /dG/ or $line =~ /ENERGY/ or eof ) {

            # process the current structure
            foreach my $p ( @pair ) {
    #            print STDERR qq{$p->[0]\t$p->[1]\n};
            }
            push @stemlist, basepairToStem(@pair);

                         #print STDERR "$stemlist[0]\n";
            my $dnb = basepairToDotsAndBrackets( \@pair, $length );
            push my @comment_list, ($dnb);
            $self->comment( \@comment_list );

            my ( $length, $dg, $equal, $energy, $id, $doc ) = split " ", $line,
              6;
            #last if ( $energy > $mfe + $ddG );    # done if energy limit reached

            $seq = uc $seq;
            $self->sequence($seq) unless ($self->sequence);    # sequence for all structures should be the same
            $seq  = "";
            @pair = ();
        }
        else {

            # add to current structure.  Initially structures are stored as an
            # array of base-paired positions
            my ( $left, $base, $pos, $x1, $right, $x2 ) = split " ", $line;

            # not sure if the following test is necessary
            # next if ( (!(defined($base))) && (!(defined($pos) )) );

            $seq .= $base;

            # in CT file, pairs are listed from both directions
            next unless ( $right && $left < $right );
            push @pair, [ $left, $right ];
        }
    }

    # remove duplicate stems and split if stems contain more than max_unpaired
    # bases

    # create unique stemlist
    my %unique;
    foreach my $stem (@stemlist) {
        my $key = $stem->left1 . "." . $stem->right1 . ".";
        $key .= $stem->left2 . "." . $stem->right2 . ".";
        $unique{$key} = $stem;
    }

    #print Dumper(\%unique ), "\n";

    my @ustem;
     my %ukey;
    foreach my $stem ( keys %unique ) {
        my @split = splitStems( $max_unpaired, $unique{$stem} );
        foreach my $ss (@split) {

            # stems must be at least 2 bp
            next if ( length( $ss->vienna_right ) < 2 );

            my $key = $ss->left1 . "." . $ss->right1 . ".";
            $key .= $ss->vienna_right . $ss->vienna_left;
            unless ( defined $ukey{$key} ) {
                $ukey{$key} = 1;
                push @ustem, $ss;
            }
        }
    }

    $self->stemList( \@ustem );

    return $self->nStem;
}

# End of mfoldMultiple
           #------------------------------------------------------------------------------
# bpseqToTopology
#
# Read the .bpseq file and return a stem list and a dots-and-bracket file
# with sequence
#
# SAMPLE FILE:
#
# Filename: a.I1.b.Bacteriophage.SP01.A2.DP.g31.bpseq
# Organism: Bacillus phage SPO1
# Accession Number: M37686
# Citation and related information available at http://www.rna.ccbb.utexas.edu
# 1 C 0
# 2 G 0
# 3 U 0
# 4 U 0
# 5 U 0
# 6 G 0
# 7 A 0
# 8 G 0
# 9 U 0
# 10 A 0
#
#
# USAGE
#    my $n_stem = $topology->bpseqToTopology( $filename, $max_unpaired );
#------------------------------------------------------------------------------
sub bpseqToTopology {
    my ( $self, $filename, $max_unpaired ) = @_;
    my $n_stem;

    open( my $in, "<", $filename )
      or die
"Topology::bpseqToToplogy unable to open input BPSEQ file ($filename)\n\n";

    my ( $organism, $access_number );
    my ( $seq, $head, @pair, @stemlist );
    while ( my $line = <$in> ) {
        next if ( $line =~ /#/ );    # skip comments
        next unless ( $line =~ /\w/ );    # skip blank lines
        chomp $line;
        ($organism) = $line =~ /^Organism: (.*)/ if ( $line =~ /^Organism/ );
        ($access_number) = $line =~ /^Accession Number: (.*)/
          if ( $line =~ /^Accession Number: / );
        if (eof) {
                      # process the current structure
            if ($head) {
                my $headn = length($head);
                print STDERR "   $headn leading N/n's removed from $filename\n";
            }

            #my ( $npair, $nseq ) = basepairPrune( \@pair, $seq );
        my $npair = \@pair;
        my $nseq = $seq;

            $self->sequence($nseq);
            my $length = length($nseq);
            $self->length($length);

            push @stemlist, basepairToStem( @{$npair} );
            my $dnb = basepairToDotsAndBrackets( $npair, $length );
            push my @comment_list, ($dnb);
            $self->comment( \@comment_list );
        }
        if ( $line =~ /^\d/ ) {

            # add to current structure.  Initially structures are stored as an
            # array of base-paired positions

            my ( $left, $base, $right ) = split " ", $line;
            if ( $base eq 'n' or $base eq 'N' ) {

                # remove leading N/n's from the sequence
                $head .= $base unless ($seq);
            }
            else {
                $seq .= $base;

                my $headn = 0;
                $headn = length($head) if ($head);
                $left  = $left - $headn;
                $right = $right - $headn;

                # in BPSEQ file, pairs are listed from both directions
                next unless ( $right && $left < $right );
                push @pair, [ $left, $right ];

            }

        }
    }
            my $id = ">$filename Organism: $organism Accession Number: $access_number";
    $self->sequenceId($id);

    # remove duplicate stems and split if stems contain more than max_unpaired
    # bases

    # create unique stemlist
    my %unique;
    foreach my $stem (@stemlist) {
        my $key = $stem->left1 . "." . $stem->right1 . ".";
        $key .= $stem->left2 . "." . $stem->right2 . ".";
        $unique{$key} = $stem;
    }

    my @ustem;
    my %ukey;
    foreach my $stem ( keys %unique ) {
        my @split = splitStems( $max_unpaired, $unique{$stem} );
        foreach my $ss (@split) {

            # stems must be at least 2 bp
            next if ( length( $ss->vienna_right ) < 2 );

            my $key = $ss->left1 . "." . $ss->right1 . ".";
            $key .= $ss->vienna_right . $ss->vienna_left;
            unless ( defined $ukey{$key} ) {
                $ukey{$key} = 1;
                push @ustem, $ss;
            }
        }
    }

    $self->stemList( \@ustem );

    return $self->nStem;
}

# End of bpseqToTopology
                        #------------------------------------------------------------------------------
# ViennaToTopology
#
#  Reads in a .va file, convert it to RNA topology.
#
# USAGE
#    $nstems = $rna->ViennaToTopology( $filename, $max_unpaired );
#------------------------------------------------------------------------------
sub ViennaToTopology {
    my ( $self, $filename, $max_unpaired ) = @_;
    open ( IN, "<$filename" ) || die "can not open $filename!\n";
    while ( my $line = <IN> ) {
    next unless ( $line =~ /\(/ or $line =~ /\)/ );
    chomp $line;
    my $vienna = $line;
    my $hash = {
        vienna => $vienna,
    };
    my $xrna = new Topology( $hash );
        $self->addTopology( $xrna );
    }
    close IN;


    # create unique stemlist
    my %unique;
    my @stemlist = @{ $self->stemList };
    foreach my $stem (@stemlist) {
        my $key = $stem->left1 . "." . $stem->right1 . ".";
        $key .= $stem->left2 . "." . $stem->right2 . ".";
        $unique{$key} = $stem;
    }

    my @ustem;
    my %ukey;
    foreach my $stem ( keys %unique ) {
        my @split = splitStems( $max_unpaired, $unique{$stem} );
        foreach my $ss (@split) {

            # stems must be at least 2 bp
            next if ( length( $ss->vienna_right ) < 2 );

            my $key = $ss->left1 . "." . $ss->right1 . ".";
            $key .= $ss->vienna_right . $ss->vienna_left;
            unless ( defined $ukey{$key} ) {
                $ukey{$key} = 1;
                push @ustem, $ss;
            }
        }
    }
                                   $self->stemList( \@ustem );
    #print $self->XIOSWrite;
    return $self->nStem;
}
# End of VieenaToTopology

       #------------------------------------------------------------------------------
# HNumToTopology
#
# Read the .h-num file and return a stem list
# The stems being read in should have a h-num value lower than $h_num_threshold
#
# SAMPLE FILE:
#
# level   length  istart  jstart  h-num
# 1       5       29      43      3.4
# 1       4       42      53      3.25
# 1       4       9       63      3.25
# 1       5       51      65      3.2
#
#
# USAGE
#    my $n_stem = $topology->HnumToTopology( $filename, $max_unpaired,
#       $h_num_threshold, $h_num_precent );
#------------------------------------------------------------------------------
sub HNumToTopology {
    my ( $self, $filename, $max_unpaired, $h_num_max, $h_num_percent ) = @_;
    my $h_num_threshold;

    my $in;
    open( $in, "<", $filename )
      or die
"Topology::bpseqToToplogy unable to open input .h-num file ($filename)\n\n";

    my @h_num;
    while ( my $line = <$in> ) {
        next if ( $line =~ /level/ );    # skip comments
        next unless ( $line =~ /\w/ );   # skip blank lines
        chomp $line;
        my ( $level, $length, $istart, $jstart, $h_num ) = split " ", $line;
        push @h_num, $h_num;
    }

    @h_num = sort { $a <=> $b } @h_num;
    my $nstem = scalar @h_num;
    print STDERR "$nstem stems read!\n";
    if ($h_num_percent) {
        # keep stems with the lowest 10% h-num values
        $h_num_threshold = $h_num[ int( @h_num * $h_num_percent / 100 ) ];
        print STDERR "Keeping $h_num_percent% stems!\n";
    }
    else {
        $h_num_threshold = $h_num_max;
    }
    print STDERR "Pruning stems using h-num threshold $h_num_threshold\n";
    close $in;
         925,14        25%
         
    open( $in, "<", $filename )
      or die
"Topology::bpseqToToplogy unable to open input .h-num file ($filename)\n\n";

    while ( my $line = <$in> ) {
        next if ( $line =~ /level/ );    # skip comments
        next unless ( $line =~ /\w/ );   # skip blank lines
        chomp $line;

        my ( $level, $length, $istart, $jstart, $h_num ) = split " ", $line;
        my $left1        = $istart;
        my $left2        = $istart + $length - 1;
        my $right1       = $jstart - $length + 1;
        my $right2       = $jstart;
        my $vienna_left  = "(" x $length;
        my $vienna_right = ")" x $length;

        if ( $h_num > $h_num_threshold ) {
            print STDERR "Skip stem [$left1 $left2 $right1 $right2] with h-num value $h_num\n";
            next;
        }

        my $stem = new Stem;
        $stem->setStem(
            [ $left1, $left2, $right1, $right2, $vienna_left, $vienna_right ] );
        $self->addStem($stem);
    }

    close $in;

    # remove duplicate stems and split if stems contain more than max_unpaired
    # bases

    # create unique stemlist
    my %unique;
    my @stemlist = @{ $self->stemList };
    foreach my $stem (@stemlist) {
        my $key = $stem->left1 . "." . $stem->right1 . ".";
        $key .= $stem->left2 . "." . $stem->right2 . ".";
        $unique{$key} = $stem;
    }

    my @ustem;
    my %ukey;
    foreach my $stem ( keys %unique ) {
        my @split = splitStems( $max_unpaired, $unique{$stem} );
        foreach my $ss (@split) {

            # stems must be at least 2 bp
            next if ( length( $ss->vienna_right ) < 2 );

            my $key = $ss->left1 . "." . $ss->right1 . ".";
            $key .= $ss->vienna_right . $ss->vienna_left;
            unless ( defined $ukey{$key} ) {
                $ukey{$key} = 1;
                push @ustem, $ss;
            }
        }
    }

    $self->stemList( \@ustem );

    return $self->nStem;
}

# End of HNumToTopology

#------------------------------------------------------------------------------
# EnergyToTopology
#
# Read the .plot energy file and return a stem list
# The stems being read in should have an energy value lower than $dG_threshold
#
# SAMPLE FILE:
#
# level   length  i       j       energy
# 1       2       321     326     -1688
# 1       2       317     325     -1679
# 1       2       315     326     -1689
# 1       2       316     323     -1689
#
# USAGE
#    my $n_stem = $topology->EnergyToTopology( $filename, $max_unpaired, $dG_threshold );
#------------------------------------------------------------------------------
sub EnergyToTopology {
    my ( $self, $filename, $max_unpaired, $dG_threshold ) = @_;
    my %dG;

    my $in;

#open( $in, "<", $filename ) or
#   die "Topology::bpseqToToplogy unable to open input .h-num file ($filename)\n\n";

    #my @dG;
    #while ( my $line = <$in> ) {
    #    next if ( $line =~ /level/ );       # skip comments
    #    next unless ( $line =~ /\w/ );  # skip blank lines
    #    chomp $line;
    #    my ( $level, $length, $istart, $jstart, $dG ) = split " ", $line;
    #    push @dG, $dG;
    #}

    #@dG = sort { $a <=> $b } @dG;
    #unless ( $dG_threshold ) {
    #    $dG_threshold = $h_num[int(@h_num * 0.65)];
    #    print STDERR "Pruning stems using delta-G threshold $dG_threshold\n"
    #}
    #close $in;

    open( $in, "<", $filename )
      or die
"Topology::bpseqToToplogy unable to open input .h-num file ($filename)\n\n";

    while ( my $line = <$in> ) {
        next if ( $line =~ /level/ );    # skip comments
        next unless ( $line =~ /\w/ );   # skip blank lines
        chomp $line;

        my ( $level, $length, $istart, $jstart, $dG ) = split " ", $line;
        my $left1        = $istart;
        my $left2        = $istart + $length - 1;
        my $right1       = $jstart - $length + 1;
        my $right2       = $jstart;
        my $vienna_left  = "(" x $length;
        my $vienna_right = ")" x $length;

        my $id = $left1 . '.' . $left2 . '.' . $right1 . '.' . $right2;
        $dG{$id} = $dG;

#if ( $dG > $dG_threshold ) {
#    print STDERR "Skip stem [ $left1 $left2 $right1 $right2 ] with delta-G value $dG\n";
#    next;
#}

        my $stem = new Stem;
        $stem->setStem(
            [ $left1, $left2, $right1, $right2, $vienna_left, $vienna_right ] );
        $self->addStem($stem);
    }

    close $in;

    # remove duplicate stems and split if stems contain more than max_unpaired
    # bases

    # create unique stemlist
    my %unique;
    my @stemlist = @{ $self->stemList };
    foreach my $stem (@stemlist) {
        my $key = $stem->left1 . "." . $stem->right1 . ".";
        $key .= $stem->left2 . "." . $stem->right2 . ".";
        $unique{$key} = $stem;
    }

            unless ( defined $ukey{$key} ) {
                $ukey{$key} = 1;
                push @ustem, $ss;
            }
        }
    }

    $self->stemList( \@ustem );

    #return $self->nStem;
    my $nstem = $self->nStem;

    return ( $nstem, \%dG );

}

# End of EnergyToTopology

#-----------------------------------------------------------------------------
# pairs
#
# read in a structure from an array of paired bases.  basically the same as
# except the paired bases are provided in a data structure of
#
#   [ [base1, base1'], [base2, base2'] ...]
#
# USAGE
#   $ n_stems  = $topology->pairs( $filename );
#-----------------------------------------------------------------------------
sub pairs {
    my ( $self, $pairs, $max_unpaired ) = @_;

    my @stemlist = basepairToStem( @{$pairs} );
    @stemlist = splitStems( $max_unpaired, @stemlist );

    my $n_stem = 0;
    foreach my $stem (@stemlist) {
        $self->addStem($stem);

        #print "stem: $n_stem\n";
        #$stem->dump;
        $n_stem++;
    }

    return ($n_stem);
}
#-----------------------------------------------------------------------------
# pairArray
#
# read in a structure from an the paired array format.  this allows an abstract
# topology to be read.  A paired array might be something like
# @p = ( 1, 5, 2, 4, 3, 6 ), indicating three stems: stem 1 pairing position 1
# and 5, stem 2 pairing postions 2 and 4, and stem 3 pairing positions 3 and 6.
#
# USAGE
#   $ n_stems  = $topology->pairArray( $ref_pair_array );
#-----------------------------------------------------------------------------
sub pairArray {
    my ( $self, $parray ) = @_;

    # stem = { 'left1'  => left start,
    #          'left2'  => left end,
    #          'right1' => right start,
    #          'right2' => right end
    #          'vienna_left' => vienna string for left side
    #          'vienna_right => vienna string for right side
    #         }
    my $n_stem = 0;

    my $pos = 0;
    while ( $pos < @$parray ) {
        my $stem = Stem->new(
            {
                left1        => $$parray[$pos],
                left2        => $$parray[$pos],
                right1       => $$parray[ $pos + 1 ],
                right2       => $$parray[ $pos + 1 ],
                vienna_left  => '(',
                vienna_right => ')'
            }
        );
        $pos += 2;
        $self->addStem($stem);
        $n_stem++;
    }

    return ($n_stem);
}

#----------------------------------------------------------------------------
# energyplot
#
# open energy dot plot file output by unafold(mfold) and read out stems, sorts
# them by energy. Stems that are in conflict (both substems) with lower-energy
# stems are dropped.
#
# USAGE
#      $n_stems = $topology->energyplot( $file );
#---------------------------------------------------------------------------
sub energyplot {
    my ( $self, $file ) = @_;

    my @stems      = ();
    my @uniq_stems = ();
    open( PLOT, "$file" ) or die "can't open plotfile $file: $!\n";
    while (<PLOT>) {
        if ( $_ =~ /\d+\s+(\d+)\s+(\d+)\s+(\d+)\s+(.+)/ ) {
            if ( $1 > 1 ) {    # skip one base stems
                    # stem: left start, right end, stem length and energy
                push @stems, [ $2, $3, $1, $4 ];
            }
        }
    }

    # sort stems by energy, left start, right end
    my @sorted_stems =
      sort { $a->[3] <=> $b->[3] || $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] }
      @stems;
    my $best_e = $sorted_stems[0][3];

    for my $newstem (@sorted_stems) {
        my $ls_new = ${$newstem}[0];
        my $le_new = ${$newstem}[0] + ${$newstem}[2] - 1;
        my $re_new = ${$newstem}[1];
        my $rs_new = ${$newstem}[1] - ${$newstem}[2] + 1;

        if ( ${$newstem}[3] == $best_e ) {
            push @uniq_stems, [ $ls_new, $le_new, $rs_new, $re_new, "", "" ];
            next;
        }
        else {
            my $flag = 0;    # flag for whether the new stem is unique or not
            for my $uniqstem (@uniq_stems) {
                if ( $ls_new <= ${$uniqstem}[1] && $le_new >= ${$uniqstem}[0] )
                {

                    # left substems overlap
                    if (   $rs_new <= ${$uniqstem}[3]
                        && $re_new >= ${$uniqstem}[2] )
                    {

                        # right substems overlap
                        $flag = 1;
                        last;
                    }
                }
            }
            if ( $flag != 1 ) { # if there is no overlap with prev stems, add it
                push @uniq_stems,
                  [ $ls_new, $le_new, $rs_new, $re_new, "", "" ];
                $flag = 0;
            }
        }
    }
    my $n_stem = 0;
    for my $stem_coor (@uniq_stems) {
        my $stem = Stem->new($stem_coor);
        $self->addStem($stem);

        #      print "@$stem_coor\n";
        #      print "stem obj: $stem\n";
        #      $stem->dump;
        #      print " stem $n_stem =>  @$stem_coor \n";
        $n_stem++;
    }

    #  print "Num stems: $n_stem\n";
    return ($n_stem);
}

#-----------------------------------------------------------------------------
# rnaml
#
# read a structure in rnaml format. Will not corrrectly read multiple
# structures in one file.
#
# USAGE
#   $n_stems = $topology->rnaml( $rnaml_file );
#-----------------------------------------------------------------------------
sub rnaml {
    my ( $self, $filename, $max_unpaired ) = @_;

    open( RNAML, "<$filename" )
      || die "Unable to open RNAML file $filename\n\n";

    my $in_sequence  = 0;
    my $in_data      = 0;
    my $in_structure = 0;

    my $sequence = "";
    my $seqlen   = 0;
    my @pairs;

    while ( my $line = <RNAML> ) {
        next if ( $line =~ /^\s*$/ );
        chomp $line;

        if ( $line =~ /<\s*rna/ ) {
            my ($title) = $line =~ /<rna\s+name\s*=\s*\"([^"]*)\"/;
            my ( $id, $doc ) = split " ", $title, 2;
            $self->sequenceId($id);
            $self->sequenceDoc($doc);
        }
        if ( $line =~ /<\s*sequence/ ) { $in_sequence = 1; }
        if ( $line =~ /<\s*data/ ) {
            if ($in_sequence) { $in_data = 1; next; }
        }
        if ( $line =~ /<\s*structure/ ) { $in_structure = 1; next; }

        # closing tags

        # read only the first structure
        if ( $line =~ /<\s*\/\s*rnaml/ )     { last; }
        if ( $line =~ /<\s*\/\s*sequence/ )  { $in_sequence = 0; }
        if ( $line =~ /<\s*\/\s*data/ )      { $in_data = 0; }
        if ( $line =~ /<\s*\/\s*structure/ ) { $in_structure = 0; }

        # data element (raw string)

        if ($in_data) {

            $line =~ s/\W//g;   # remove non-alphabetic characters from sequence
            $sequence .= $line;
            $seqlen = length($sequence);

  # structure elements, e.g.,  <base-pair _5p-base-num="2"  _3p-base-num="403"/>

        }
        elsif ($in_structure) {

            my ( $l, $r ) = $line =~
/base-pair\s*_5p-base-num\s*=\s*\"(\d+)\"\s*_3p-base-num\s*=\s*\"(\d+)\"/;
            unless ( defined($l) && defined($r) ) {
                print "undefined l or r:$line\n";
            }

            if ( $l < $r ) {
                push @pairs, [ $l, $r ];
            }

        }
    }
    close RNAML;

    $self->sequence($sequence);
    my @stemlist = basepairToStem(@pairs);
    @stemlist = splitStems( $max_unpaired, @stemlist );

    my $n_stem = 0;
    foreach my $stem (@stemlist) {
        $self->addStem($stem);

        #print "stem: $n_stem\n";
        #$stem->dump;
        $n_stem++;
    }
    $self->sequence($sequence);

    return ($n_stem);
}

# end of rnaml


"""
