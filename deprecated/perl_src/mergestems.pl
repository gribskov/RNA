#!/apps/group/bioinformatics/apps/perl-5.16.1/bin/perl -w
#!/usr/bin/perl -w 

#------------------------------------------------------------------------------
# $Id: mergestems.pl,v 2.1 2014/07/22 15:29:41 gribskov Exp $
# 
# Merge topologically similar stems.  All of the merge cases in this version 
# corresponde to the newMergeCases of Reazur's plotstem_newrules
#
# Copied from Jiajie's version of plotstem_newrules 4 Apr 2012
#
# Revision log at end of file
#------------------------------------------------------------------------------
#use strict;
use List::Util qw(min max);
use Getopt::Std;
use GD;
use lib 'RNA';  # use a symbolic link to set this to your prefered perl_src
use lib "/group/mgribsko/rna/RNA/src/perl_src/";
use lib "/apps/group/bioinformatics/apps/perl-5.16.1/bin/perl/";
use Topology;

my $USAGE = qq{mergestem.pl [-c <mergecases>] [-g <ddG>] [l <int>] [-p <plot_file>]  [-q] [-m <curated_file>] <CT_file> };

my $STEM_CENTER_THRESHOLD = 6;
my $DEF_DELTA_DELTAG      = 5;

my %option;
getopts( 'c:g:hl:p:m:q', \%option );

# help
if ( $option{h} ) {
    print STDERR "$USAGE\n\n";
    exit 1;
}

# which mergecases should be run, and in what order
#my $casenum = '123333345123451234512345';
#my $casenum = '12345';
my $casenum = '';


if ( $option{c} ) {
    $casenum = $option{c};
} 
my $curated_file;
if ( $option{m} ) {
    $curated_file = $option{m};
}

# delta deltaG
my $ddG = $DEF_DELTA_DELTAG;
if ( $option{g} ) {
    $ddG = $option{g};
}

# specify sequence length (may not be the same as the real sequence length)
my $seqlen = 400;
if ( $option{l} ) {
    $seqlen = $option{l};
}

# output file for plots, no plot if undefined
my $plotfile = '';
if ( $option{p} ) {
    $plotfile = $option{p};
    # make sure the file name ends in"." so 1.png and 2.png can be added
    unless ( $plotfile =~ /\.$/ ) {
        $plotfile =~ s/\.*\s*$/./;
    }
}

# quiet - do not write messages to stderr - better for web version
my $quiet = 0;
if ( $option{q} ) {
    $quiet = 1;
}

# input file
my $filename = $ARGV[0];
unless ( -e $filename ) {
    print STDERR "   Unable to open input CT file ($filename)\n\n";
    exit 2;
}

#print STDERR "curated_file: $curated_file ||| input file: $filename\n"; 

#die;

# read stems

my $rna = new Topology;
my $nstem = $rna->mfoldMultiple( $filename, 3, $ddG );
print STDERR "   $nstem stems read from $filename\n" unless $quiet;;

unless ( $option{l} ) {
    # do not override -l option
    $seqlen = length( $rna->sequence ) + 1;
}
my $stemlist = $rna->stemListAsArray;
my $cstems = stemCenters( $stemlist );
my $ncent = @{$cstems};

#my $jum = 0;
#foreach my $st (@{$cstems}) {
#    print "st num: $jum; @$st\n";
#    $jum++;
#}

# tree contraction
# initialize tree with a dummy node that includes the whole sequence as the
# left coord, right coord and stem center.  coordinates begin at 1 and end at
# sequence length  

my @nodelist;           # list of all nodes in tree
my $tree = treeNode( 0, [ 1, 1, $seqlen, $seqlen, ($seqlen+1)/2 ] ); 
$nodelist[0] = $tree;

# add remaining nodes to tree
my $ntree = 0;
foreach my $s ( @{$cstems} )  {
    $ntree++;
    my $node = treeNode( $ntree, $s );
    $nodelist[$ntree] = $node;
    treeAddNode( $node, $tree );
}
my $nnode = @nodelist;
treeAddParent( $tree );
#treeDump( $tree );

if ( $plotfile ) {    
    # plot original nodes as PNG    
    my @plotlist = @nodelist;
    shift @plotlist;        # remove root node for plotting
    my $plot     = plotInit( $seqlen, $ntree );
    my $nplotted = plotStems( $plot, \@plotlist );
    my $file1    = $plotfile."1.png";
    plotPrint( $plot, $file1 );
}

# the merging is specified by the string $casenum.  merging is performed in the 
# order indicated by the string, e.g., the string "152" indicates run three 
# merge cases, in the order 1, 5, 2.

$casenum =~ s/[^\d]//g;     # remove non digits
my @cases = split "", $casenum;

my $nmerge = 0;
while ( @cases ) {
    my $case = shift @cases;
    $nmerge += mergeCase1(    @nodelist ) if ( $case == 1 );
    $nmerge += mergeCase2(    @nodelist ) if ( $case == 2 );
    $nmerge += mergeCase3( 4, @nodelist ) if ( $case == 3 );
    #larger threshold does a better job of merging multiple nested stems for mergecase 3
    #$nmerge += mergeCase3( 8, @nodelist ) if ( $case == 3 );
    $nmerge += mergeCase4( 4, @nodelist ) if ( $case == 4 );
    $nmerge += mergeCase5( 8, @nodelist ) if ( $case == 5 );

}
print STDERR "\n   $nmerge stems merged\n" unless $quiet;

# save the merged tree to a topology for printing.  After merging, any stem with
# an 'x' in its ID has been merged

my $merged = new Topology;
my $plotlist;
for my $i ( 1 .. $#nodelist ) {
#    print "ID: $nodelist[$i]->{id}\n";
#    print "$nodelist[$i]->{pos}->[0]    $nodelist[$i]->{pos}->[1]     $nodelist[$i]->{pos}->[2]    $nodelist[$i]->{pos}->[3]\n";
#    if (($nodelist[$i]->{pos}->[0] == 278) && ($nodelist[$i]->{pos}->[1] == 286)) { #	print "XXXXXXXXX\n";  }
    if ( $nodelist[$i]->{id} !~ /x/ ) {
        my $s = new Stem;
        $s->left1 ( $nodelist[$i]->{pos}->[0] );
        $s->left2 ( $nodelist[$i]->{pos}->[1] );
        $s->right1( $nodelist[$i]->{pos}->[2] );
        $s->right2( $nodelist[$i]->{pos}->[3] );
        $merged->addStem( $s );

        push @{$plotlist}, $nodelist[$i];       
    } 
}

$nnode = @{$plotlist};
print STDERR "   $nnode nodes after pruning\n\n" unless $quiet;
#treeDump( $tree );


if (( $plotfile ) && ($curated_file) ) { 
    my $ctopology_obj = new Topology;
    my $num_curated_stem = $ctopology_obj->parseManualCuratedFile($curated_file);
    
    my $cur_stemlist = $ctopology_obj->stemListAsArray;
    my $num = 0;    
    foreach my $stem ( @{$cur_stemlist} ) {
	$num++;
	my $center = 0;
        foreach my $point ( @{$stem} ) {
            $center += $point;
        }
	$center = $center / 4.0;
	my $new_cstem = [@{$stem}, $center ];
	my $id = $num.'c'; #The curated stems have 'c', concatenated to their id
#	print STDERR "curated stem added: $id: @$new_cstem\n";
	my $node = treeNode( $id, $new_cstem );
	push @{$plotlist}, $node; 
    } # all curated stems have been added
    @{$plotlist} = sort { $b->{pos}->[3] - $b->{pos}->[0] <=> $a->{pos}->[3] - $a->{pos}->[0] || 
			      $a->{pos}->[0] <=> $b->{pos}->[0] } @{$plotlist}; #sort by stem width and start position
    
    # print out plot as PNG, curated stems are added to the plot
    my $nstem    = @{$plotlist};
    my $plot     = plotInit( $seqlen, $nstem );
    my $nplotted = plotStemsWithCurated( $plot, $plotlist );
    my $file2    = $plotfile."2.png";
    plotPrint( $plot, $file2 );
    
} elsif ( $plotfile ) {
    
    # print out plot as PNG
    my $nstem    = @{$plotlist};
    my $plot     = plotInit( $seqlen, $nstem );
    my $nplotted = plotStems( $plot, $plotlist );
    my $file2    = $plotfile."2.png";
    plotPrint( $plot, $file2 );
}

my $xiosfile = $filename.".xios";
#my $merged->XIOSWrite($xiosfile);
print $merged->XIOSWrite;

exit 0;

#------------------------------------------------------------------------------
# stemFilter
#
# removes stems that are exactly identical from the list.  Should not be needed 
# any more.  this is done in mfoldMultiple.
#
# USAGE
#    $filtered_stems = stemFilter( $stems );
#------------------------------------------------------------------------------
sub stemFilter{
    my ( $stems ) = @_;
    
    my $MINLEN = 1;
    
    my $newstems;
    push @{$newstems}, [ @{$stems->[0]} ];
    
    foreach my $s ( @{$stems} ) {
       
       my $same  = 0;
       foreach my $n ( @{$newstems} ) {
           
           my $different = 0;
           foreach my $i ( 0 .. $#{$stems->[0]} ) {
               if ( $s->[$i] != $n->[$i] ) {
                   $different = 1;
                   last;
               }
           }
           unless ( $different ) {
               $same = 1;
               last;
           }
       }
       unless ( $same ) {
           if ( $s->[1] - $s->[0] > $MINLEN and $s->[3] - $s->[2] > $MINLEN ) {
               push @{$newstems}, [ @{$s} ];
           }
       }    
    }
    
    return $newstems;
}

# End of stemFilter

#------------------------------------------------------------------------------
# stemCenters
#
# adds the location of the stem center to the stem object.  the stem center is 
# calculated as the average of the beginning and ending positions of both left 
# and right half-stems.  This is different than the original diffgraph method 
# which used only the loop edges.
#
# USAGE
#    $sorted_centered_stems = stemCenters( $stems );
#------------------------------------------------------------------------------
sub stemCenters{
    my ( $stems ) = @_;    
    my $cstems;
    
    foreach my $stem ( @{$stems} ) {
        my $center = 0;
        foreach my $point ( @{$stem} ) {
            $center += $point;
        }
        $center = $center / 4.0;
        push @{$cstems}, [ @{$stem}, $center ];
    }
    
    $cstems = [ sort { $b->[3]-$b->[0] <=> $a->[3]-$a->[0] ||
                       $a->[0]         <=> $b->[0]            } 
                @{$cstems} 
              ];
    
    return $cstems;
}

# End of stemCenters

#------------------------------------------------------------------------------
# plotInit
#
# sets up a canvas for plotting stem positions.
#
# USAGE
#    $plot = plotInit( $seqlen );
#------------------------------------------------------------------------------
sub plotInit{
    my ( $seqlen, $nstem ) = @_;
    
    my $plot = {};
    my $hcanvas = 1550;
    my $vcanvas = 1200;
    
    $plot->{hcanvas} = $hcanvas;
    $plot->{vcanvas} = $vcanvas;
    $plot->{vpad}    = 100;
    $plot->{hpad}    = 100;
    $plot->{hscale}  = ( $hcanvas - 2*$plot->{hpad} ) / ($nstem+2);
    #print "seqlen = $seqlen\n";
    $plot->{vscale}  = ( $vcanvas - 2*$plot->{vpad} ) / $seqlen;
    
    my $im = new GD::Image( $hcanvas, $vcanvas) || die;
    my $white = $im->colorAllocate( 255, 255, 255 );
    my $black = $im->colorAllocate( 0, 0, 0 );
    my $red   = $im->colorAllocate( 255, 0, 0 );
    my $blue  = $im->colorAllocate( 0, 0, 255 );
    my $green = $im->colorAllocate(0, 100, 0); #dark green, for light green (0, 255, 0)

    my $dash = $im->setStyle( 1,1,1,1,0,0,);
    $im->setThickness( 1.0 );
    #my $arial = "/usr/openwin/lib/locale/ru.ansi-1251/X11/fonts/TrueType/ArialRegular.ttf";
    
    # basic box and gridlines
            
    $im->rectangle( $plot->{hpad}, 
                    $plot->{vpad}, 
                    $plot->{hcanvas}-$plot->{hpad}, 
                    $plot->{vcanvas}-$plot->{vpad}, 
                    1 );
                    
    my $pos = 10;
    while ( $pos < $seqlen ) {
        
        $im->setStyle( 1,1,1,1,1,1,0,0,);
        if ( $pos % 50 ) {
           $im->setStyle( 1,0,0,0,0,0) 
        }
        $im->line( $plot->{hpad},
                   $plot->{vcanvas} - $plot->{vpad} - $pos*$plot->{vscale},
                   $plot->{hcanvas} - $plot->{hpad},
                   $plot->{vcanvas} - $plot->{vpad} - $pos*$plot->{vscale},
                   "gdStyled" );
        $pos += 10;    
    }
    $plot->{im} = $im;
        
        
    return $plot;
}

# End of plotInit

#------------------------------------------------------------------------------
# plotStems
#
# plot the stems.  If there are merged stems, they are shown as offset boxes on 
# the parent (merged) stem.
#
# USAGE
#    plotStems( );
#------------------------------------------------------------------------------
sub plotStems{
    my ( $plot, $stems ) = @_;
    my $OVERLAP = 3;
    my $nstem = 1;
    my $center_old = -100;
    my $smallstep = $plot->{hscale}/10.0;
    
    foreach my $stem ( @{$stems} ) {
        my ( $left1, $left2, $right1, $right2, $center ) = @{$stem->{pos}};
        my $x1 = $nstem * $plot->{hscale} + $plot->{hpad};
        my $x2 = ( $nstem+1) * $plot->{hscale} + $plot->{hpad};
        my $y1 = $plot->{vcanvas} - $plot->{vpad} - $left1 * $plot->{vscale};
        my $y2 = $plot->{vcanvas} - $plot->{vpad} - $left2 * $plot->{vscale};
        $plot->{im}->filledRectangle( $x1, $y2, $x2, $y1, 2 ); #both half stem are colored green
        $plot->{im}->rectangle( $x1, $y2, $x2, $y1, 1 ); 
        $plot->{im}->string( "gdSmallFont", $x1, $y1, $stem->{id}, 1 ); #print id
        $plot->{im}->string( "gdSmallFont", $x1, $y2, $nstem, 1 ); # merge stem num ber, depends on how the sorting is done
	
        my $y3 = $plot->{vcanvas} - $plot->{vpad} - $right1 * $plot->{vscale};
        my $y4 = $plot->{vcanvas} - $plot->{vpad} - $right2 * $plot->{vscale};
        $plot->{im}->filledRectangle( $x1, $y4, $x2, $y3, 3 );
        $plot->{im}->rectangle( $x1, $y4, $x2, $y3, 1 );
    
        # merged stems 
        if ( @{$stem->{merge}} > 1 ) {
            foreach my $mstem ( @{$stem->{merge}} ) {
		my ( $id, $left1, $left2, $right1, $right2, $center ) = @{$mstem};
                $x1 += $smallstep;
                $x2 = $x1 + $plot->{hscale};
		
                my $y1 = $plot->{vcanvas} - $plot->{vpad} - $left1 * $plot->{vscale};
                my $y2 = $plot->{vcanvas} - $plot->{vpad} - $left2 * $plot->{vscale};
                $plot->{im}->rectangle( $x1, $y2, $x2, $y1, 1 ); 
                
                my $y3 = $plot->{vcanvas} - $plot->{vpad} - $right1 * $plot->{vscale};
                my $y4 = $plot->{vcanvas} - $plot->{vpad} - $right2 * $plot->{vscale};
                $plot->{im}->rectangle( $x1, $y4, $x2, $y3, 1 );
	    }
        }
	# line connecting boxes
	$x1 = $nstem * $plot->{hscale} + $plot->{hpad};
        $x1 += $plot->{hscale}/2.0;
        $plot->{im}->line( $x1, $y3, $x1, $y2, 1 );
	
        $nstem++; 
    }
    
    return $nstem;
}

# End of plotStems

#------------------------------------------------------------------------------
# plotStemsWithCurated
#
# plot the predicted stem along with curated.  If there are merged stems, they are shown as offset boxes on 
# the parent (merged) stem.
#
# USAGE
#    plotStemsWithCurated( );
#------------------------------------------------------------------------------
sub plotStemsWithCurated {
    my ( $plot, $stems ) = @_;
    my $OVERLAP = 3;
    my $nstem = 1;
    my $center_old = -100;
    my $smallstep = $plot->{hscale}/10.0;
    
    foreach my $stem ( @{$stems} ) {
        my ( $left1, $left2, $right1, $right2, $center ) = @{$stem->{pos}};
	if (  $stem->{id} =~ /c/ ) {  #if the stem is curated stem
	    my $x1 = $nstem * $plot->{hscale} + $plot->{hpad};
	    my $x2 = ( $nstem+1) * $plot->{hscale} + $plot->{hpad};
	    my $y1 = $plot->{vcanvas} - $plot->{vpad} - $left1 * $plot->{vscale};
	    my $y2 = $plot->{vcanvas} - $plot->{vpad} - $left2 * $plot->{vscale};
	    $plot->{im}->filledRectangle( $x1, $y2, $x2, $y1, 4 ); #both half stem are colored green
	    $plot->{im}->rectangle( $x1, $y2, $x2, $y1, 1 ); 
	    $plot->{im}->string( gdSmallFont, $x1, $y1, $stem->{id}, 1 ); #print id
	    #$plot->{im}->string( gdSmallFont, $x1, $y2, $nstem, 1 ); # merge stem num ber, depends on how the sorting is done
	    
	    my $y3 = $plot->{vcanvas} - $plot->{vpad} - $right1 * $plot->{vscale};
	    my $y4 = $plot->{vcanvas} - $plot->{vpad} - $right2 * $plot->{vscale};
	    $plot->{im}->filledRectangle( $x1, $y4, $x2, $y3, 4 ); #both half stem are colored green
	    $plot->{im}->rectangle( $x1, $y4, $x2, $y3, 1 );
	    # line connecting boxes
	    $x1 = $nstem * $plot->{hscale} + $plot->{hpad};
	    $x1 += $plot->{hscale}/2.0;
	    $plot->{im}->line( $x1, $y3, $x1, $y2, 1 );
	    
	} else { # for all the other predicted stems
	    
	    my $x1 = $nstem * $plot->{hscale} + $plot->{hpad};
	    my $x2 = ( $nstem+1) * $plot->{hscale} + $plot->{hpad};
	    my $y1 = $plot->{vcanvas} - $plot->{vpad} - $left1 * $plot->{vscale};
	    my $y2 = $plot->{vcanvas} - $plot->{vpad} - $left2 * $plot->{vscale};
	    $plot->{im}->filledRectangle( $x1, $y2, $x2, $y1, 2 );
	    $plot->{im}->rectangle( $x1, $y2, $x2, $y1, 1 ); 
	    $plot->{im}->string( gdSmallFont, $x1, $y1, $stem->{id}, 1 ); #print id

	    my $y3 = $plot->{vcanvas} - $plot->{vpad} - $right1 * $plot->{vscale};
	    my $y4 = $plot->{vcanvas} - $plot->{vpad} - $right2 * $plot->{vscale};
	    $plot->{im}->filledRectangle( $x1, $y4, $x2, $y3, 3 );
	    $plot->{im}->rectangle( $x1, $y4, $x2, $y3, 1 );

	    $plot->{im}->string( gdSmallFont, $x1, $y3, $nstem, 1 ); # show merge stem numberin figure, depends on how the sorting is done
	    	    
	    #------------------------------------------------------------

	    # merged stems 
	    if ( @{$stem->{merge}} > 1 ) {
		foreach my $mstem ( @{$stem->{merge}} ) {
		    
		    my ( $id, $left1, $left2, $right1, $right2, $center ) = @{$mstem};
		    $x1 += $smallstep;
		    $x2 = $x1 + $plot->{hscale};
		    
		    my $y1 = $plot->{vcanvas} - $plot->{vpad} - $left1 * $plot->{vscale};
		    my $y2 = $plot->{vcanvas} - $plot->{vpad} - $left2 * $plot->{vscale};
		    $plot->{im}->rectangle( $x1, $y2, $x2, $y1, 1 ); 
		    
		    my $y3 = $plot->{vcanvas} - $plot->{vpad} - $right1 * $plot->{vscale};
		    my $y4 = $plot->{vcanvas} - $plot->{vpad} - $right2 * $plot->{vscale};
		    $plot->{im}->rectangle( $x1, $y4, $x2, $y3, 1 );
		    
		}
	    }
	    
	    # line connecting boxes
	    
	    $x1 = $nstem * $plot->{hscale} + $plot->{hpad};
	    $x1 += $plot->{hscale}/2.0;
	    $plot->{im}->line( $x1, $y3, $x1, $y2, 1 );
	    
	} #end else
	
        $nstem++; 
    } #end foreach my $stem ( @{$stems} ) 
    
    return $nstem;
}

# End of plotStems

#------------------------------------------------------------------------------
# plotPrint
#
# print the plot to the output file.
#
# USAGE
#    plotPrint( $plot, $plotfile );
#------------------------------------------------------------------------------
sub plotPrint{
    my ( $plot, $plotfile ) = @_;
        
    open( PLT, ">$plotfile" ) or die "Unable to open plot file ($plotfile)\n";
    binmode PLT;
    print PLT $plot->{im}->png(8);
    close PLT;
        
    return;
}

# End of plotPrint

################################################################################
# The merge cases detect various redundant relationships between a parent, a, 
# and child stem, b and, when detected, merge the child stem into the parent 
# stem.  The parent stem has coordinates (a0, a1, a2, a3) and the child stem 
# (b0, b1, b2, b3)
################################################################################

#------------------------------------------------------------------------------
# mergeCase1
#
#   |aaaaaaaaaa|        |aaaaaaaaaa|
#     |bbbbbb|            |bbbbbb|
#
# merge children that have left and right half stems that are both contained in 
# the halfstems of their parent
#
#   a0 <= b0 && b1<=a1 && a2 <= b2 && b3< a3
#
# USAGE
#   $nmerged = mergeCase1( $nodelist );
#------------------------------------------------------------------------------
sub mergeCase1 {    
    my ( @nodelist ) = @_;
    my $nmerged = 0;
    
    # stem center test on a threshold
    # checks to see if merging the child and parent would violate it's 
    # relationship with other stems

    my @to_be_merged = ();    
    for my $i ( @nodelist ) {
        # skip already merged stems and root node
        next if ( $i->{id} =~ /x/ || $i->{id} == 0 );   
        
        my ($a0, $a1, $a2, $a3, $ac) = ( $i->{pos}->[0], $i->{pos}->[1], 
                                         $i->{pos}->[2], $i->{pos}->[3], 
                                         $i->{pos}->[4]
                                       );
        
        foreach my $child ( @{$i->{children}} ) {
            next if ( $child->{id} =~ /x/ );    # skip alread merged stems
            
            my ($b0, $b1, $b2, $b3, $bc) = ( $child->{pos}->[0], $child->{pos}->[1], 
                                             $child->{pos}->[2], $child->{pos}->[3], 
                                             $child->{pos}->[4]
                                           );
        
            if ( $a0 <= $b0 && $b1 <= $a1 && $a2 <= $b2 && $b3 <= $a3 ) {
                
                printf STDERR "   Case 1: Child node %s [ %d %d %d %d ] merged into ", 
                              $child->{id},$b0, $b1, $b2, $b3 unless $quiet;
                
                printf STDERR "%s [ %d %d %d %d ]\n",
                              $i->{id},$a0, $a1, $a2, $a3, $ac  unless $quiet;                
                push @to_be_merged, [ $i, $child ]; # child is merged to parent
                $nmerged++;
            }
        }
    }
    
    # merge all the nodes, not sure why this isn't done inside the loop as for
    # other cases
    foreach my $node ( @to_be_merged ) {
        treeMergeNodes( @$node );
    }        
    
    return $nmerged;
}

# end of mergeCase1

#------------------------------------------------------------------------------
# mergeCase2
#
#   |0aaaaaaaa1|      |2aaaaaaaa3|
#         |0bbbb1|   |2bbbb3|
#
# merges when 
# 1) the left half stem of the parent partially/completely overlaps the left 
#    halfstem of the child, AND  
# 2) the right half stem of the parent partially/completely overlaps left 
#    halfstem of the child, AND
# 3) the stem center difference should be less than $STEM_CENTER_THRESHOLD
#
# in this situation, the child has some bases that lie outside the parent half 
# stems on the right hand side.  if the child halfstems lie entirely within the
# parent they will be merged by an earlier invocation of mergeCase1.
#
# a is the parent stem and b is the child
#
# ( $b0 >= $a0 && $b0 <= $a1 && $b1 < $a2 )     &&
# ( $b3 >= $a2 && $b3 <= $a3 && $b2>$a1 )       && 
# ( $b1 > $a1 || $b2 < $a2 )                    &&
# ( $stemcenter_diff <= $STEM_CENTER_THRESHOLD )
# 
# USAGE
#   $nmerged = mergeCase2( $nodelist );
#   Note: uses global variable $STEM_CENTER_THRESHOLD
#------------------------------------------------------------------------------
sub mergeCase2{
    my ( @nodelist ) = @_;
    
    my $nmerged = 0;
    for my $i ( @nodelist ) {
        # skip already merged stems and the root node
        next if ( $i->{id} =~ /x/ || $i->{id} == 0 );
        
        # the parent stem
        my ( $a0, $a1, $a2, $a3, $ac ) = ( $i->{pos}->[0], $i->{pos}->[1], 
                                           $i->{pos}->[2], $i->{pos}->[3], 
                                           $i->{pos}->[4]
                                         );
        
        my @children = @{$i->{children}};
        foreach my $child ( @children ) {
            next if ( $child->{id} =~ /x/ );    # skip already merged
            # the child stem
            my ( $b0, $b1, $b2, $b3, $bc ) = ( $child->{pos}->[0], $child->{pos}->[1], 
                                               $child->{pos}->[2], $child->{pos}->[3], 
                                               $child->{pos}->[4]
                                             );
            
            if ( ($b0 >= $a0 && $b0 <= $a1 && $b1 < $a2) &&
                 ($b3>= $a2 && $b3<=$a3 && $b2>$a1)      && 
                 ($b1 > $a1 || $b2 < $a2)                   ) {
                # condition for case 2, has been met, 
                # Is there any need to test for unmatched bases? I think not...
                # TODO: check that the loop region is large enough to be a real loop?
                
                my $stemcenter_diff = abs($bc - $ac);
                if ( $stemcenter_diff <= $STEM_CENTER_THRESHOLD ) {
                    printf STDERR "   Case 2:Child node %s [ %d %d %d %d ] merged into ", 
                    $child->{id},$b0, $b1, $b2, $b3 unless $quiet;
                    
                    printf STDERR "%s [ %d %d %d %d ]\n",
                    $i->{id},$a0, $a1, $a2, $a3, $ac  unless $quiet;
                    
                    treeMergeNodes( $i, $child );
                    $nmerged++;
                }
            }
        }
    }
    return $nmerged;
}

# end of mergeCase2

#------------------------------------------------------------------------------
# mergeCase3
#
#   |0aaaaaaaa1|                    |2aaaaaaaa3|
#                |0bbbb1|  |2bbbb3|
#
# merges child and parent nodes when 
# 1) the child is completely nested within the parent stem, AND
# 2) the gap between child and parent stems is not greater than $gap_limit, AND
# 3) the centers are less than $STEM_CENTER_THRESHOLD apart
#
# ( $b0 > $a1 && $b1 < $a2 && $b2 > $a1 && $b3 < $a2 ) and 
#                 $left_halfstem_gap <= $threshold and 
#                 $right_halfstem_gap <= $threshold and
#                 $stemcenter_diff <= $STEM_CENTER_THRESHOLD
#
# ( $b0 > $a1 && $b1 < $a2 && $b2 > $a1 && $b3 < $a2 ) &&                
# ( $b0-$a1 <= $threshold && $a2-$b3 <= $threshold )   &&
# ( $abs($bc - $ac) <= $STEM_CENTER_THRESHOLD ) 
# 
# $b0 - $a1 is the gap between the parent and child left half stems
# $a2 - $b3 is the gap between the parent and child right half stems
# abs($bc - $ac) is the difference in the center position of the two stems
#
# USAGE
#   $nmerged = mergeCase3( $gap_limit, $nodelist );
#   Note: uses global variable $STEM_CENTER_THRESHOLD
#------------------------------------------------------------------------------
sub mergeCase3{
    my ( $threshold, @nodelist ) = @_;
    
    my $nmerged = 0;
    for my $i ( @nodelist ) {
        # skip already merged stems and the root node
        next if ( $i->{id} =~ /x/ || $i->{id} == 0 );
        
        # the parent stem
        my ( $a0, $a1, $a2, $a3, $ac ) = ( $i->{pos}->[0], $i->{pos}->[1], 
                                           $i->{pos}->[2], $i->{pos}->[3], 
                                           $i->{pos}->[4]
                                         );
        
        my @children = @{$i->{children}};
        foreach my $child ( @children ) {
            next if ( $child->{id} =~ /x/ );
            my ( $b0, $b1, $b2, $b3, $bc ) = ( $child->{pos}->[0], $child->{pos}->[1], 
                                               $child->{pos}->[2], $child->{pos}->[3], 
                                               $child->{pos}->[4]
                                             );
            
            if ( $b0 > $a1 && $b1 < $a2 && $b2 > $a1 && $b3 < $a2 ) {
                # don't need to check that a loop is physically reasonable since
                # the nested stem has a loop
                
                my $stemcenter_diff    = abs($bc - $ac);
                my $left_halfstem_gap  = $b0 - $a1;
                my $right_halfstem_gap = $a2 - $b3;
                
                if ( $left_halfstem_gap  <= $threshold and 
                     $right_halfstem_gap <= $threshold and
                     $stemcenter_diff    <= $STEM_CENTER_THRESHOLD ) {
                    
                    printf STDERR "   Case 3: Child node %s [ %d %d %d %d ] merged into ", 
                    $child->{id}, $b0, $b1, $b2, $b3 unless $quiet;
                    
                    printf STDERR "%s [ %d %d %d %d ]\n",
                    $i->{id},$a0, $a1, $a2, $a3, $ac  unless $quiet;
                    
                    treeMergeNodes( $i, $child );
                    $nmerged++;
                }
            }
        }
    }
    return $nmerged;
}

# end of mergeCase3

#------------------------------------------------------------------------------
# mergeCase4
#
#       |0aaaaaaaa1|           |2aaaaaaaa3|
#        |0bbbb1|     |2bbbb3|
#
#       |0aaaaaaaa1|           |2aaaaaaaa3|
#                    |0bbbb1|     |2bbbb3|
#
# merges nested child node with parent nodes when 
# 1) one child halfstem overlaps partially/completely with the parent AND
# 2) the other half stem is nested inside the parent half stems AND
# 3) the nested half stem is less than $gap_limit away from the boundaries of 
#    the parent AND
# 4) the stem center difference is less than $STEM_CENTER_THRESHOLD
#
# ( $b0 <= $a1 && $b0 >= $a0 && $b3 <  $a2 && $b2 >  $a1 ) || 
# ( $b0 > $a1  && $b1 < $a2  && $b3 >= $a2 && $b3 <= $a3 )
# Either the left halfstem overlap or the right halfstem 
#
# USAGE
#   $nmerged = mergeCase3( $gap_limit, $nodelist );
#   Note: uses global variable $STEM_CENTER_THRESHOLD
#------------------------------------------------------------------------------
sub mergeCase4{
    my ( $threshold, @nodelist ) = @_;
    
    my $nmerged = 0;
    for my $i ( @nodelist ) {
        # skip already merged stems and the root node
        next if ( $i->{id} =~ /x/ || $i->{id} == 0 ); # what does this mean? $i->{id} == 0
        # the parent stem
        my ( $a0, $a1, $a2, $a3, $ac ) = ( $i->{pos}->[0], $i->{pos}->[1], 
                                           $i->{pos}->[2], $i->{pos}->[3], 
                                           $i->{pos}->[4]
                                         );
        
        my @children = @{$i->{children}};
        foreach my $child ( @children ) {
            next if ( $child->{id} =~ /x/ );
            my ( $b0, $b1, $b2, $b3, $bc ) = ( $child->{pos}->[0], $child->{pos}->[1], 
                                               $child->{pos}->[2], $child->{pos}->[3], 
                                               $child->{pos}->[4]
                                             );
            
            # The left halfstem overlaps, but the right does not, implicit b1< a2
            if ($b0 <= $a1 && $b0 >= $a0 && $b3 < $a2 && $b2 > $a1 ) {
                my $stemcenter_diff = abs($bc - $ac);
                my $right_halfstem_gap = $a2 - $b3;
                if ( $right_halfstem_gap <= $threshold and
                     $stemcenter_diff <= $STEM_CENTER_THRESHOLD ) {
                    
                    printf STDERR "   Case 4: Child node %s [ %d %d %d %d ] merged into ", 
                    $child->{id}, $b0, $b1, $b2, $b3 unless $quiet;
                    
                    printf STDERR "%s [ %d %d %d %d ]\n",
                    $i->{id},$a0, $a1, $a2, $a3, $ac  unless $quiet;
                    
                    treeMergeNodes( $i, $child );
                    $nmerged++;
                }
            }
            
            # The right halfstem overlaps, but the left does not, implicit b2> a1
            if ($b0 > $a1 && $b1 < $a2 && $b3 >= $a2 && $b3 <= $a3) {
                my $stemcenter_diff = abs($bc - $ac);
                my $left_halfstem_gap = $b0 - $a1;
                if ( $left_halfstem_gap <= $threshold and
                     $stemcenter_diff <= $STEM_CENTER_THRESHOLD ) {
                    
                    printf STDERR "   Case 4: Child node %s [ %d %d %d %d ] merged into ", 
                    $child->{id}, $b0, $b1, $b2, $b3 unless $quiet;
                    
                    printf STDERR "%s [ %d %d %d %d ]\n",
                    $i->{id}, $a0, $a1, $a2, $a3, $ac  unless $quiet;
                    
                    treeMergeNodes( $i, $child );
                    $nmerged++;
                }
            }
        }
    }
    return $nmerged;
}

# end of mergeCase4

#------------------------------------------------------------------------------
# mergeCase5
#
#       |0aaaaaaaa1|                    |2aaaaaaaa3|
#     |0bbbb1|                            |2bbbb3|
#
#       |0aaaaaaaa1|                    |2aaaaaaaa3|
#          |0bbbb1|                           |2bbbb3|

#
# Looks for stems that overlap and are both nested within the same parent.  For
# this case a and b are children of the same parent, d.
#
# merge non-nested stems, when both halfstem overlaps
# number of matched and unmatched bases are below a threshold 
#
# The stem center difference also has to be below a threshold
# 
# ( $b0 < $a0 && $b1 >= $a0 && $b3 <= $a3 && $b3 >= $a2 ) || 
# ( $b3 > $a3 && $b2 <= $a3 && $b0 >= $a0 && $b0 <= $a1 )     
# plus, satisfy, matched, unmatched and stem center threshold
#
# USAGE
#   $nmerged = mergeCase5( $threshold, $nodelist );
#------------------------------------------------------------------------------
sub mergeCase5{
    my ( $threshold, @nodelist ) = @_;
    
    my $nmerged = 0;
    for my $i ( @nodelist ) {
        # skip already merged stems
        next if ( $i->{id} =~ /x/ );
        
        # parent stem
        my ( $d0, $d1, $d2, $d3, $dc ) = ( $i->{pos}->[0], $i->{pos}->[1], 
                                           $i->{pos}->[2], $i->{pos}->[3], 
                                           $i->{pos}->[4]
                                         );
        
        my @children = @{$i->{children}};
        foreach my $c_one ( 0 .. $#children ) {
            my $c1 = $children[$c_one];
            next if ( $c1->{id} =~ /x/ );   # skip already merged stems
            
            # child 1 stem
            my ( $a0, $a1, $a2, $a3, $ac ) = ( $c1->{pos}->[0], $c1->{pos}->[1], 
                                               $c1->{pos}->[2], $c1->{pos}->[3], 
                                               $c1->{pos}->[4]
                                             );
            
            my ($left_unmatched, $left_matched)   = (0, 0);
            my ($right_unmatched, $right_matched) = (0, 0);
            
            foreach my $c_two ( $c_one+1 .. $#children ) {
                my $c2 = $children[$c_two];
                next if ( $c2->{id} =~ /x/ );   # skip already merged stems
                
                # child 2 stem
                my ( $b0, $b1, $b2, $b3, $bc ) = ( $c2->{pos}->[0], $c2->{pos}->[1], 
                                                   $c2->{pos}->[2], $c2->{pos}->[3], 
                                                   $c2->{pos}->[4]
                                                 );
                
                #
                # 
                if ( ( $b0 < $a0 && $b1 >= $a0 && $b3 <= $a3 && $b3 >= $a2 ) || 
                     ( $b3 > $a3 && $b2 <= $a3 && $b0 >= $a0 && $b0 <= $a1 )    ) { 
                    
                    # _in and _out are inner and outer boundries of the halfstem
                    # there are two possible cases depending on whether a is to 
                    # the left or to the right of b
                    
                    my $l1_in = $b0; my $l1_out = $a0;                  # left halfstem, left inner and outer boundries of the halfstem
                    if ( $a0 > $b0 ) { $l1_in = $a0; $l1_out = $b0; }   
                    
                    my $l2_in = $b1; my $l2_out = $a1;                  # left halfstem, right inner and outer boundries of the halfstem
                    if ( $a1 < $b1 ) { $l2_in = $a1; $l2_out = $b1; }

                    my $r1_in = $b2; my $r1_out = $a2;                  # right halfstem, left inner and outer boundries of the halfstem
                    if ( $a2 > $b2 ) { $r1_in = $a2; $r1_out = $b2; }

                    my $r2_in = $b3; my $r2_out = $a3;                  # right halfstem, right inner and outer boundries of the halfstem
                    if ( $a3 < $b3 ) { $r2_in = $a3; $r2_out = $b3; }

                    # different length stems would create a imbalance in the
                    # unmatched base count, so it needs to be normalized                  
                    my $left_stem_len_diff  = abs(($a1-$a0) - ($b1-$b0));
                    my $right_stem_len_diff = abs(($a3-$a2) - ($b3-$b2));
                    $left_unmatched   = $l1_in - $l1_out + $l2_out - $l2_in; 
                    $left_unmatched  -= $left_stem_len_diff;
                    
                    $right_unmatched  = $r1_in - $r1_out + $r2_out - $r2_in,
                    $right_unmatched -= $right_stem_len_diff;
                    
                    my $c1_len = int(($a1-$a0 + 1 + $a3 - $a2 + 1)/2); #average of stem length
                    my $c2_len = int(($b1-$b0 + 1 + $b3 - $b2 + 1)/2);
                    my $min_stem_len = 0;
                    $min_stem_len = $c2_len;
                    if ( $c1_len < $c2_len ) { $min_stem_len = $c1_len; } #smaller of the two average of stem length
                                       
                    $left_matched  = $l2_in - $l1_in + 1; #atleast two bases have to be matched
                    $right_matched = $r2_in - $r1_in + 1;
                    # it can also determined from stem length or 30% (which ever 
                    # is larger). the unmatched parameter become a bit tricky, 
                    # when merging a long stem with a smaller stem.
                    my $MIN_MATCH_THRESHOLD = 1; # for stem length of 3,2
                    if ( $min_stem_len > 3 ) {
                        $MIN_MATCH_THRESHOLD = 2; # it remains 2, for 4,5,6
                        if ( 0.3 * $min_stem_len > 2 ) { 
                            $MIN_MATCH_THRESHOLD = int( 0.3*$min_stem_len ); # 1/3 of average stem length 
                        }
                    } 
                    
                    # Unmatched threshold can be 6
                    # unmatch threshold can be 50% of half stem length or which 
                    # ever is smaller... Times 2 or 8
                    $threshold = $min_stem_len; # 50% of base pairs can be unmatched
 
                    
                    my $stemcenter_diff = abs($bc - $ac);
                    
                    if ( $left_unmatched  <= $threshold           and  
                         $left_matched    >= $MIN_MATCH_THRESHOLD and
                         $right_unmatched <= $threshold           and 
                         $right_matched   >= $MIN_MATCH_THRESHOLD and
                         $stemcenter_diff <= $STEM_CENTER_THRESHOLD + 1 ) {
                        
                        printf STDERR "   Case 5: Child node %s [ %d %d %d %d ]merged to",
                                      $c1->{id}, $a0, $a1, $a2, $a3, $ac  unless $quiet;
                        printf STDERR "child node %s [ %d %d %d %d ]\n",
                                      $c2->{id}, $b0, $b1, $b2, $b3 unless $quiet;
                        
                        treeMergeNodes( $c1, $c2 );
                        $nmerged++;
                    }                
                }
            }
        }
    }
    return $nmerged;
}

# end of mergeCase5

#------------------------------------------------------------------------------
# treeNode
#
# defines a node in the stem tree structure.  Each node is an hashref.
#
# USAGE
#    $new_tree_node = tree_node( $stem );
#------------------------------------------------------------------------------
sub treeNode{
    my ( $id, $stem ) = @_;
    my $node = {};
    
    $node->{id} = $id;
    $node->{pos} = [ @{$stem} ];
    $node->{children} = [];
    $node->{neighbor} = [];
    $node->{done} = 0;
    $node->{maxgap} = 0;

    # push the node onto merge list so that its original coordinates are 
    # preserved
    
    push @{$node->{merge}}, [ $id, @{$stem} ];     
    
    return $node;
}
        
# end of treeNode

#------------------------------------------------------------------------------
# treeAddNode
#
# since stems are presented in order of their loopsize it is impossible for a 
# later stem to be completely outside an earlier one - at most it can overlap.
# this makes one less case to deal with.
#
# USAGE
#   treeAddNode( $node, $tree );
#------------------------------------------------------------------------------
sub treeAddNode{
    my ( $node, $tree ) = @_;
    
    my @stack;
        
#    printf "new node %s [%d %d %d %d %d]\n", 
#           $node->{id},
#           $node->{pos}->[0],
#           $node->{pos}->[1],
#           $node->{pos}->[2],
#           $node->{pos}->[3],
#           $node->{pos}->[4] unless $quiet;
        
    push @stack, $tree;
    while ( @stack ) {
        my $tnode = shift @stack;
        # check each of the children, if this node is included in one, follow each
        #  branch.  otherwise add this node as a sibling

        my $sibling = 1;
        foreach my $child ( @{$tnode->{children}} ) {
            if ( treeInside( $node, $child ) ) {
                $sibling = 0;
                
                unless ( $child->{done} ) {
                    push @stack, $child;
                    $child->{done} = 1;
                }
                        
            }
        }
                
        # $sibling is true if the new node was not inside any of the children
        # in this case add it to the child list of the current node
                
        if ( $sibling ) {
            # check to make sure that this node has not already been added
            # as a child of this node (remember it is a reticulated tree)
            
            #TODO check that this test is not needed
            
#            my @childlist = @{$tnode->{children}};
#            unless ( $childlist[$#childlist] == $node ) {
        #    if ( @{$tnode->{children}} ) {
        #        unless ( $tnode->{children}->[$#{$tnode->{children}}] == $node ) {
                    
        #                   push @{$tnode->{children}}, $node;
    #                       printf STDERR "   adding as child [%d %d %d %d]\n", 
    #                                 $node->{pos}->[0],
    #                                 $node->{pos}->[1],
    #                                 $node->{pos}->[2],
    #                                 $node->{pos}->[3];
        #        } else {
                    #print STDERR "   node already added\n";
        #        }
        #    } else {
                push @{$tnode->{children}}, $node;
        #    }
                }      
        }
        treeResetDone( $tree );
                        
        return $tree;
}

# end of treeAddNode

#------------------------------------------------------------------------------
# treeAddParent
#
# Add a back-pointer to the parent(s) of each node.
#
# USAGE
#    treeAddParent( $tree );
#------------------------------------------------------------------------------
sub treeAddParent{
        my ( $tree ) = @_;
        
        my @stack;
        push @stack, $tree;

        
        # traverse tree making a back link from each child to its parent
        
        while ( @stack ) {
            my $node = pop @stack;
            
            unless ( $node->{done} ) {
            foreach my $child ( @{$node->{children}} ) {                
                push @{$child->{parent}}, $node;
                push @stack, $child;
            }
            }
        $node->{done} = 1;
        }
        
        
        treeResetDone( $tree );
        
        return;
}

# End of treeAddParent

#------------------------------------------------------------------------------
# treeInside
#
# returns true if the first node is completely inside the second.  That is if
# left1-1 >left2-2 and right2-1 < right1-2
# where the -x suffix indicates node 1 or node 2
#
# USAGE
#    $is_inside = treeInside( $node1, $node2 );
#------------------------------------------------------------------------------
sub treeInside{
    my ( $node1, $node2) = @_;
    my $inside = 0;
        
    my $stemcenter_diff = abs($node1->{pos}->[4] - $node2->{pos}->[4]);
        # test if node1 is inside boundary of node 2, overlap is ok if the 
        # centers are within $threshold distance
        if ( $node1->{pos}->[0] >= $node2->{pos}->[0] and       
             $node1->{pos}->[3] <= $node2->{pos}->[3] and 
             $stemcenter_diff <= $STEM_CENTER_THRESHOLD) {
            # check if first node(node1) spills over into second node(node2), 
            # i.e. node1 left-halfstem extends into node2 right-halfstem, and 
            # vice versa 
        if ( $node2->{pos}->[1] < $node1->{pos}->[2] && 
             $node2->{pos}->[2] > $node1->{pos}->[1]    ) {
             $inside = 1;
        }
    }
        
    return $inside;
}

# End of treeInside

#------------------------------------------------------------------------------
# treeSameChildren
#
# checks to see if the children lists are the same with the exception of the
# nodes themselves (i.e., they are the same if the only difference is b is 
# nested in a).
#
# USAGE
#    $is_the_same = treeSameChildren( $node1, $node2 );
#------------------------------------------------------------------------------
sub treeSameChildren{
    my ( $node1, $node2 ) = @_;
    
    # if neither node has children they don't match    
    unless ( @{$node1->{children}} and @{$node2->{children}} ) {
        return 0;
    }
    
    # make a hash of the union of children of each node
    my $is_same = 1;    
    my %list;
    foreach my $child ( @{$node1->{children}} ) {
        unless ( $child != $node1 and $child->{id} =~ /x/ ) {
           $list{$child} = 1;
        }
    }
    
    foreach my $child ( @{$node2->{children}} ) {
        unless ( $child != $node2 and $child->{id} =~ /x/ ) {
           $list{$child}++;
        }
    }

    # if there are no children, we can quit
    unless ( keys %list ) {
        return 0;
    }
        
    # check to see all the children occur twice
    foreach my $child ( keys %list ) {
        if ( $list{$child} < 2 ) {
            # children that only occur once are not in common.  if there is even 
            # one the lists are not the same
            $is_same = 0;
            last;
        }
    }
        
    return $is_same;
}

# End of treeSameChildren

#------------------------------------------------------------------------------
# treeParentMatch
#
# checks the parents of the two specified nodes.  If one parent set is the 
# of the other, the node with the superset is returned as parent and the node
# with the subset is returned as the child.
#
# USAGE
#    ( $parent, $child ) = treeParentMatch( $node1, $node2 );
#------------------------------------------------------------------------------
sub treeParentMatch{
    my ( $node1, $node2 ) = @_;
    
    # make the node with more parents the parent and the other the child
    
    my %list;
    my $parent = $node1;
    my $child = $node2;
    if ( @{$node2->{parent}} > @{$node1->{parent}} ) {
        $parent = $node2;
        $child = $node1;
    }
        
    # hash of parents bigger (parent) node, omitting child if present
    
    foreach my $p ( @{$parent->{parent}} ) {
        next if ( $p->{id} =~ /x/ );
        unless ( $p eq $child ) {
            $list{$p}++;
        }
    }
    
    # check parents of child (smaller) node, if any are not defined, 
    # child is not a subset
    
    foreach my $p ( @{$child->{parent}} ) {
        next if ( $p eq $parent or $p->{id} =~ /x/ );
        unless ( defined $list{$p} ) {
            $parent = undef;
            $child = undef;
            last;
        }
    }
    
    return ( $parent, $child );
}

# End of treeParentMatch

#------------------------------------------------------------------------------
# treeGetParent
#
# return a list of non-merged parents
#
# USAGE
#    @parentlist = treeGetParent( $node );
#------------------------------------------------------------------------------
sub treeGetParent{
    my ( $node ) = @_;
    my @parents;
    
    foreach my $p ( @{$node->{parent}} ) {
        unless ( $p->{id} =~ /x/ ) {
            push @parents, $p;
        }
    }
    
    return @parents;
}

# End of treeGetParent

#------------------------------------------------------------------------------
# treeGetChild
#
# return a list of non-merged children of this node
#
# USAGE
#    @childlist = treeGetChild( $node );
#------------------------------------------------------------------------------
sub treeGetChild{
    my ( $node ) = @_;
    my @childlist;
    
    foreach my $c ( @{$node->{children}} ) {
        unless ( $c->{id} =~ /x/ ) {
            push @childlist, $c;
        }
    }
    
    
    return @childlist;
}

# End of treeGetChild

#------------------------------------------------------------------------------
# treeParentCount
#
# Counts the number of non-merged parents
#
# USAGE
#    $count = treeParentCount( $node );
#------------------------------------------------------------------------------
sub treeParentCount{
    my ( $node ) = @_;
    my $count = 0;
    
    foreach my $p ( @{$node->{parent}} ) {
        unless ( $p->{id} =~ /x/ ) {
            $count++;
        }
    }
    
    return $count ;
}

# End of treeParentCount

#------------------------------------------------------------------------------
# treeChildCount
#
# return the number of non-merged children
#
# USAGE
#    $count = treeChildCount( $node );
#------------------------------------------------------------------------------
sub treeChildCount{
    my ( $node ) = @_;
    my $count = 0;
    
    foreach my $c ( @{$node->{children}} ) {
        unless ( $c->{id} =~ /x/ ) {
            $count++;
        }
    }
    
    return $count;
}

# End of treeChildCount

#------------------------------------------------------------------------------
# treeMergeNodes
#
# Merges the source node into the destination node.  The coordinates of the stem
# in $dnode->{pos} are updated, and all child and merge node information of the 
# source node are copied into the child and merged node information of the 
# destination node
#
# USAGE
#   $nnodes = treeMergeNodes( $dest_node, $source_node );
#------------------------------------------------------------------------------
sub treeMergeNodes{
    my ( $dnode, $snode ) = @_;
    
    $snode->{id} .= 'x';   # mark source node as merged
    
    # merge snode into dnode
    #push @{$dnode->{merge}}, [ $snode->{id}, @{$snode->{pos}} ];
    
    # copy all of the previously merged nodes in snode into the dnode merged
    # node list
    foreach my $merged_node ( @{$snode->{merge}} ) {
        push @{$dnode->{merge}}, $merged_node;
    }
    
    # make a hash of children of the destination node
    my %d_child_list;
    foreach my $child ( @{$dnode->{children}} ) {
        $d_child_list{$child}++;
    }
    
    # copy all children of snode into dnode
    foreach my $child ( @{$snode->{children}} ) {
        unless ( defined $d_child_list{$child} ) {
            # only add if it is not already a child of dnode
            push @{$dnode->{children}}, $child;
        }
        
        # if dnode is not already a parent of child, add it
        my $is_parent = 0;
        foreach my $cp ( @{$child->{parent}}) {
            if ( $cp eq $dnode ) { $is_parent = 1; }
        }
        unless ( $is_parent ) {             
            push @{$child->{parent}}, $dnode; 
        }
    }
    
    # update dnode pos
    
    if ( $snode->{pos}->[0] < $dnode->{pos}->[0] ) { $dnode->{pos}->[0] = $snode->{pos}->[0]; }
    if ( $snode->{pos}->[1] > $dnode->{pos}->[1] ) { $dnode->{pos}->[1] = $snode->{pos}->[1]; }
    if ( $snode->{pos}->[2] < $dnode->{pos}->[2] ) { $dnode->{pos}->[2] = $snode->{pos}->[2]; }
    if ( $snode->{pos}->[3] > $dnode->{pos}->[3] ) { $dnode->{pos}->[3] = $snode->{pos}->[3]; }
    my $nnodes = @{$dnode->{merge}};
    
    return $nnodes;
}

# end of treeMergeNodes

#------------------------------------------------------------------------------
# treeNodeList
#
# Return an array of nodes that can be traced as children of the starting node
#
# USAGE
#    $nodelist = treeNodeList( treenode );
#------------------------------------------------------------------------------
sub treeNodeList{
    my ( $tree ) = @_;
    
    my $nodelist = [];
    
    my @stack;
    push @stack, $tree;
    treeResetDone( $tree );
    
    # traverse tree adding each child to list
    
    while ( @stack ) {
        my $node = pop @stack;
        next if ( $node->{done} );
        
        push @{$nodelist}, $node;       
        push @stack, @{$node->{children}};      
        $node->{done} = 1;
    }
    treeResetDone( $tree );
    
    @{$nodelist} = sort { $a->{pos}->[0]<=>$b->{pos}->[0] } @{$nodelist};
    
    return $nodelist;
}

# End of treeNodeList

#------------------------------------------------------------------------------
# treeResetDone
#
# marks all nodes in the subtree of tree with $node->{done} = 0.  This resets
# the done flag for all nodes in the subtree and makes ready for new processing.
#
# USAGE
#    treeResetDone( );
#------------------------------------------------------------------------------
sub treeResetDone{
    my ( $tree ) = @_;
    
    my @stack;
    push @stack, $tree;
    
    while ( @stack ) {
        my $node = pop @stack;
        $node->{done} = 0;
        
        foreach my $child ( @{$node->{children}} ) {
            next unless ( $child->{done} );
            push @stack, $child;
        }
    }
    
    return;
}

# End of treeResetDone

#------------------------------------------------------------------------------
# treeDump
#
# Prints a formatted version of the entire tree (or subtree).  Do not use this
# to generate the output of the merging; use Topology.pm instead.  This is 
# intended only for use in debugging.
#
# USAGE
#   treeDump( $tree );
#------------------------------------------------------------------------------
sub treeDump{
    my ( $tree ) = @_;
    
    my @stack;
    
    if ( $tree->{done} ) {
        treeResetDone( $tree );
    }
    $tree->{done} = 1;

    print STDERR "root node\n   children:  ";
    if ( @{$tree->{children}} ) {
        foreach my $c (0 .. $#{$tree->{children}} ) {
            my $child = $tree->{children}->[$#{$tree->{children}} - $c];
            print STDERR "$child->{id} ";
            push @stack, $child;
        }
    }
    print "\n\n";
    
    my $node = pop @stack;
    
    while ( $node ) {
        print STDERR "node:$node  id:$node->{id}  done:$node->{done}\n";
        if ( $node->{done} ) {
            $node = pop @stack;
            next;
        }
        
        printf "\nnode %2s [ %d %d %d %d ]\n",
               $node->{id}, $node->{pos}->[0], $node->{pos}->[1],
               $node->{pos}->[2], $node->{pos}->[3]; 
            
        if ( @{$node->{merge}} > 1 ) {
            printf STDERR "        [ %d %d %d %d ]\n", 
                   $node->{merge}->[0]->[1], $node->{merge}->[0]->[2],
                   $node->{merge}->[0]->[3], $node->{merge}->[0]->[4];      
            
            print STDERR "   merged: ";
            foreach my $m ( @{$node->{merge}} ) {
                printf STDERR "%s ", $m->[0] unless $quiet;
            }
             print STDERR "\n";
        }

         print STDERR "   parent: ";
        foreach my $p ( @{$node->{parent}} ) {
             print STDERR "$p->{id} ";
        }
         print STDERR "\n";
    
        if ( @{$node->{children}} ) {
            print STDERR "   children: ";
            foreach my $c (0 .. $#{$node->{children}} ) {
                my $child = $node->{children}->[$#{$node->{children}} - $c];
                print STDERR "$child->{id}\t";
                unless ( $child->{done} ) {
                        push @stack, $child;
                }
            } 
            print STDERR "\n"; 
        }
    
        print STDERR "   done = ", $node->{done},"\n\n";
        $node->{done} = 1;       
        $node = pop @stack;   
    }
    print STDERR "\n";
        
    return;
}

# End of treeDump

#-------------------------------------------------------------------------------
# treeStemRelationship
#
# return a letter indicating the relationship between node1 and node2 from the
# point of view of node 1.  For instance, if node2 is within node1, I is 
# returned
#
# usage
#	$relationship = treeStemRelationship(  node1, node2 );
#-------------------------------------------------------------------------------
sub treeStemRelationship{
	my ($node1, $node2 ) = @_;
	my $relationship = "";
	
	return $relationship;
}

# end of treeStemRelationship

#------------------------------------------------------------------------------
# $Log: mergestems.pl,v $
# Revision 2.1  2014/07/22 15:29:41  gribskov
# changed library path to a link.  fixed error in command line option syntax
# (-h has no argument)
#
# Revision 2.0  2014/07/08 15:14:46  gribskov
# Updated vs mrg121003 branch, renamed current versions to be 2.0
#
# Revision 1.15  2014/07/08 12:43:15  gribskov
# Merged with v 1.11.2.2.
# CVg Committing in .
#
# Revision 1.14  2012/12/21 07:26:42  rahmanr
#
# 1. change to merge case 5: added +1 to to both left and right matched bases to get the correct number
# 2. if a curated file is provided with a -m flag, then it plotted along with the merged stem files
#
# Revision 1.13  2012/10/03 11:55:19  gribskov
# Rolled back version after merging
#
# Revision 1.11.2.2  2012/09/21 14:20:37  gribskov
# Changed unpaired residue limit for mfoldMultiple to 3, which means that
# bulges/bubbles of 1-2 bases are allowed.
#
# Revision 1.11.2.1  2012/07/15 22:16:17  gribskov
# fixed stemCenters and treeDump which were broken by Jiajie
# added initial shell for calculating revised relationship tree, but will not proceed before
# checking whether the fix to stemcenter corrects merging problems
#
# Revision 1.11  2012/07/05 19:13:18  huang147
# changed the order of stems, now the stems with similar centers are put together. makes the plot easier to understand.
#
# Revision 1.12  2012/07/15 22:55:31  gribskov
# Reverting to version 1.9.  Jiajie has apparently been using this as her
# personal copy.
#
# Revision 1.9  2012/06/01 16:37:09  gribskov
# Modified for more complete merging
# 1) changed default scheme to 123333345 (from 12345)
# 2) changed STEM_CENTER_THRESHOLD from 3 to 6
# 3) changed gap threshold for mergecase3 from 3 to 8
# 4) fixed automatic setting of length for plot
# 5) fixed bug that deleted all data for plot 2
#
# Revision 1.8  2012/04/25 04:35:34  rahmanr
# added unless $quiet to one of the print statement
#
# Revision 1.7  2012/04/05 20:17:06  gribskov
# Modified printing to format more neatly in fingerprint_pipeline.pl
#
# Revision 1.6  2012/04/05 17:40:10  gribskov
# Small changes to printing, removed trace of nodes as they are added.
#
# Revision 1.5  2012/04/05 17:34:08  gribskov
# Updated documentation and formatting throughout.  No functional changes
# except to treeDump where all the printing code was commented out so that
# the function did nothing at all.  Commented out the call to treeDump in the
# main program so that nothing would print.  This would have been the right
# way to do it to start with.
#
# Revision 1.4  2012/04/05 13:03:42  gribskov
# updating documentation and reformatting merge cases.  Not done yet.  No
# functional changes
#
# Revision 1.3  2012/04/05 11:48:00  gribskov
# Added print of number of stems loaded.  Updated treeNode documentation and
# reformatted treeAddNode.
#
# Revision 1.2  2012/04/05 00:36:26  gribskov
# Seems to basically work.  Need to check that the mergecases are applied in
# the correct order and fix the feedback messages.
#
# Revision 1.1  2012/04/05 00:06:15  gribskov
# Initial but nearly working version.  merged stems are not yet stored back into
# a topology for printing. Still needs much more cleanup.
#
#------------------------------------------------------------------------------
