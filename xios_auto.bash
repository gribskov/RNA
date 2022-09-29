#!/bin/bash
#-------------------------------------------------------------------------------
# automatically run xios_from_rnastructure over a set of fast files with a grid
# of windows and deltaG.  deltaG range is handled by xios_from_rnastructur.py
# and only one run of Fold is required per window
#-------------------------------------------------------------------------------
family=                             # prefix on fastafile names, blank gives all
fasta=../data/fasta                 # directory with input fasta files
RNAstructure=../../RNAstructure     # directory with installed RNAstructure pkg
xioscode=../../RNA                  # directory with XIOS RNA pkg
log=autot3.log                      # log file with full interactive trace

# window range
wbeg=1
wend=15
#wend=3
winc=1

# delta deltaG range. Because bash can't do floating point math, these should be
# double the real values
dbeg=1
dend=8
#dend=2
dinc=1

for (( window=$wbeg; window<=$wend; window+=$winc )); do
   steps=1
   end=8
   dg=$(printf '%.1f,%.1f,%.1f' `echo $dbeg / 2.0 | bc -l` `echo $dend / 2.0 | bc -l` `echo $dinc / 2.0 | bc -l`)
   steps=$(( steps+1 ))
   cmd="python $xioscode/xios_from_rnastructure.py -i $fasta -r $RNAstructure -y $xioscode -f $family*.fa  -w $window -d $dg"
   echo $cmd
   $cmd &>> $log
done
