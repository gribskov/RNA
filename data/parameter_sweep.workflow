# workflow for calculating fingerprints from sequence using RNAstructure:partition program
# framework for optimizing hyperparameters, options filled by parameter_sweep.py are indicated by < > below
#       folding temperature (partition)
#       fraction of structure (stochastic_to_xios)
#       minimum bases in stem (stochastic_to_xios)
# THIS WORKFLOW WILL NOT RUN UNLESS THE VALUES IN < > ARE PROVIDED
---
  project: <$base/project_name>
  definitions:
    # paths
    python: python
    base: /scratch/bell/mgribsko/rna/
    ps: $base/work/260306_parameter_sweep
    RNAstructure: $base/RNAstructure/exe
    XIOSexe: $base/RNA
    XIOSdata: $base/RNA/data
    fasta: $XIOSexe/fasta_fixed
  directories:
    # directories to create under project; also creates symbols
    partition: $project/partition
    stochastic: $project/stochastic
    xios: $project/xios
    fpt: $project/fpt
  commands:
    partition:
      command: $RNAstructure/partition $in $out $temperature $option
      temperature: <-t 280>
      option: -q
      in: $fasta/*.fa
      out: $partition/%in.replace('.fa', '.pfs')
    stochastic:
      command: $RNAstructure/stochastic $in $out $seed
      seed: <-s 3>
      in: $partition/*.pfs
      out: $stochastic/%in.replace('.pfs', '.ct')
    xios:
      command: $python $XIOSexe/stochastic_to_xios.py $stemsize $fraction $in $out
      stemsize: <-m>
      fraction: <-c>
      in: $stochastic/*.ct
      out: $xios/%in.replace('.ct', '.xios')
    fingerprint:
      command: $python $XIOSexe/fingerprint_random.py -m $XIOSdata/2to7stem.mdb.pkl $option -r $in -f $out
      option: -l 500000 -c 3
      in: $xios/*.xios
      out: $fpt/%in.replace('.xios', 'fpt')



