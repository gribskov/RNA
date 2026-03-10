#!/bin/bash
################################################################################
# rename slurm job files (.slurm), workflow (.workflow) and stderr and stdout 
# files (.err and .out) to the results directory for each parameter set
# parameter conditions start with p_
################################################################################
for d in p_*; do
    if [[ "$d" == *\.* ]];then
        com="mv $d $thisdir"
        echo -e "\t$com"
        $com
    else
        thisdir=$d
        echo "directory $d"
    fi
done
