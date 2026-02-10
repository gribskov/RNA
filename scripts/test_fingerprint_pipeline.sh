#!/bin/bash
#-------------------------------------------------------------------------------
# this shell script tests the pipeline of programs needed to create fingerprints
# from Fasta files using the RNAstructure package. The series of programs is
# RNAstructure/partition
# RNAstructure/stochastic
# RNA/stochastic_to_xios
#
# usage:
# bash test.sh ../../standards/fasta_fixed/rnasep_a2.Alcaligenes_eutrophus.fa
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
$partition_command

# stochastic
echo -e "\nStochastic"
ct=${pfs/.pfs/.ct}
echo "ct: $ct"
stochastic_command="$rs_exe/stochastic $pfs $ct"
echo $stochastic_command
$stochastic_command

# stochastic_to_xios
echo -e "\nstochastic_to_xios"
xios=${ct/.ct/.xios}
source activate   /home/mgribsko/.conda/envs/cent7/5.3.1-py37/rna
s2rcommand="python ../../RNA/stochastic_to_xios.py $ct $xios"
echo $s2rcommand
$s2rcommand

# fingerprint_random
echo -e "\nfingerprint_random"
fpt=${xios/.xios/.fpt}
fingerprint_command="python ../../RNA/fingerprint_random.py -m ../../RNA/data/2to7stem.mdb.pkl -r $xios -f $fpt"
echo $fingerprint_command
$fingerprint_command
