#!/bin/bash
#-------------------------------------------------------------------------------
# this shell script tests the pipeline of programs needed to create fingerprints
# from Fasta files using the RNAstructure package. The series of programs is
# RNAstructure/partition
# RNAstructure/stochastic
# RNA/stochastic_to_xios
#-------------------------------------------------------------------------------
export DATAPATH=../../RNAstructure/data_tables/
rs_exe=../../RNAstructure/exe
echo "DATAPATH:$DATAPATH"
echo "RNAstructure executables:$rs_exe"

# input file
f=$*
echo "fasta: $f"

echo -e "\nPartition"
pfs=${f##*/}
pfs=${pfs/.fa/.pfs}
echo "partition function save file: $pfs"
partition_command="$rs_exe/partition $f $pfs"
echo $partition_command
#$partition_command

# stochastic
echo -e "\nStochastic"
ct=${pfs/.pfs/.ct}
echo "ct: $ct"
stochastic_command="$rs_exe/stochastic $pfs $ct"
echo $stochastic_command
#$stochastic_command

# stochastic_to_xios
echo -e "\nstochastic_to_xios"
xios=${ct/.ct/.xios}
source activate   /home/mgribsko/.conda/envs/cent7/5.3.1-py37/rna
command="python ../../RNA/stochastic_to_xios.py $ct $xios"
echo $command
$command
