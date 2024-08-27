# RNA XIOS Fingerprint

## Installation
1. Create a directory for the package
```
mkdir rna
cd rna
```
2. Install the XIOS/fingerprint package by cloning the git trepository. 
This will create a directory ***RNA*** with the XIOS/fingerprint code
```commandline
git clone https://github.com/RNA.git
```
2. Install Matthews RNAstructure package (command line version).
See https://rna.urmc.rochester.edu/RNAstructureDownload.html for more information. 
This will create a directory called ***RNAstructure***. Most of the scripts assume 
that the ***RNAstructure*** directory is a sub-directory of the parent directory of 
***RNA***, *i.e.* ***../RNAstructure***.
```
wget https://rna.urmc.rochester.edu/Releases/current/RNAstructureLinuxTextInterfaces64bit.tgz
tar -xvzf RNAstructureLinuxTextInterfaces64bit.tgz
```

## Making the motif library
Available soon

## Constructing fingerprints from sequence
* Sequences should be in FASTA format and should contain only RNA bases ACGU
* Code locations below are with respect to the ***RNA*** directory.

Multiple steps can be run using *manager_new.py* (see *Job Manager*, below), which manages 
submitting multiple jobs 
in a workflow as individual *SLURM* jobs. See *parameter_sweep.py* for an example of how 
to run *manager_new.py* on a series of three computational steps. 


### Calculate XIOS graphs from structural ensemble
#### ../RNAstructure/partition
Calculate the partition function from the FASTA file. 
The partition function is used to sample in structures according to their free energy of folding
(ΔG) in the next step.
 * input: FASTA sequence - *.fa
 * output: partition function - *.pfs
 * options
   - -t \<int> folding temperature (K), e.g. 260 - 340K

#### ../RNAstructure/stochastic
Probabalistically sample structures from the partition function according to free energy of folding.
Sampling uses a random number generator which must be seeded; the same seed will always generate the
same samples. For testing, you may want to use the same seed to guarantee that results are comparable.
The seed should be a large integer without obvious small factors - the current system 
time is often used.
 * input: *.pfs
 * output: *.ct
 * options
   - -s \<int> random number seed for sampling (default=1234 which is not a good seed)

#### stochastic_to_xios.py ####
Calculate XIOS graph from sampled ct files 
 * input: *.ct
 * output: *.xios
 * options
   - -c \<int> minimum number of counts to include paired bases as a stem (default=50)
   - -m \<int> minimum number of paired bases in a stem (default=3)

## Job manager ##
A large scale analysis, or even one with only a few hundred sequences, requires executing workflows 
of multiple programs from the RNAstructure and XIOS RNA packages. If executed as individual jobs on 
a compute cluster, thousands of jobs may need to be launched and tracked. On computers with large 
multi-processor compute nodes, the clusters we use typically use have 128-256 cpus, it is more 
efficient to run a single job that internally manages the usage of the cpus. We provide a flexible job 
manager in *new_manager/manager.py* for running arbitrary pipelines in a ***SLURM*** controlled 
cluster environment.

Features
1. Can run pipelines with any number of steps using a ***YAML*** formatted file to define the workflow
2. Allows symbols to be used to define workflow steps.
3. Wildcard inputs are used to generate separate processes for each file after filename expansion
3. All workflow steps are logged with timestamps.
3. When processes are terminated due to expiration of the time limit, the workflow can automatically
restart from the last successful step.

### Workflow file ###
The *workflow file* controls the steps in the workflow contolling inputs, outputs, and options at each step. 
The *workflow file* is divided into two main sections. ***YAML*** requires that the file begin with *---*. \
Indentation is used to group material into sections. Here's an example of a workflow with just the section headers.
```
---
# Example workflow skeleton
  definitions:
    # Add definitions here

  stage:
    # Add workflow stages/steps here
```
#### Definitions Section ####
The definitions section defines variables that can be used in the rest of the workflow. This saves a lot of typing
and makes it easier to modify workflows. All values that are likely to be changed during normal usage should
be set up as *definitions* so that the *stage* information does not need to be edited. Here's an example 
of some definitions. 
```
  definitions:
    project: human
    python: python
    fasta: ../data
    base: $project
    ct: $base/ct
    xios: ../RNA
    RNAstructure: ../RNAstructure
```
The *definiton section* above defines 7 variables that can be used in other definitions or in *stage* definitions.
Variables are referenced as *$variable* similarly to shell scripts. For instance, in the example 
above, the variable *project* is used to define the variable *base* (making *base* equal to "human").

#### Stage Section ####
The *stage section* is used to define actual commands to be executed. Each stage is entered in linear
order of execution (*i.e.* the first stage listed is executed first, the second stage listed is 
executed second, ...). This workflow defines three stages named *partition*, *stochastic*, and *xios*, but 
you could name them anything. For each starting file, *partition* is run first, then *stochastic*, 
and finally *xios*.

**command:** defines the actual command that will be executed after substituting the defined values 
for variabls ($in, $out, $option).

**in:** defines the input to the program being run. If *$in* contains wildcards, they will be 
expanded to the set of matching files and run one at a time.

**out:** defines the output file for the program. One often needs to create an output file name
from the input file name. This can be done with the *replace* operator.<br/>
***out: $partition/%in.replace('.fa', '.pfs')*** 
means that the output file should appear in the *partition/* directory and be the same as *$in*, the input file, with *.fa* replaced by *.pfs*.
```
---
  definitions:
#   code locations
    python: python
    RNAstructure: ../rna/RNAstructure
    xiosexe: ../rna/RNA

#   project input/output directories
    fasta: ../standards/fasta_fixed
    project: 240807sweep
    partition: $project/partition
    stochastic: $project/stochastic
    # program variables
    xios: $project/xios
    XTemp: 280
    XSeed: 409875209874285475029872084757
    XStem: 3
    XFraction: 50
 
  stage:
    partition:
      command: $RNAstructure/partition  $in $out $option
      option:
        -t $XTemp
        -q
      in: $fasta/*.fa
      out: :$partition/%in.replace('.fa', '.pfs')
    
    stochastic:
      command: $RNAstructure/stochastic $in $out $option
      option: -s $XSeed
      in: $partition/*.pfs
      out: :$stochastic/%in.replace('.pfs', '.ct')
    
    xios:
      command: $python $xiosexe/stochastic_to_xios.py $option $in $out
      option: -m $XStem -c $XFraction
      in: $stochastic/*.ct
      out: :$xios/%in.replace('.ct', '.xios')
```


## Libraries - file to import (classes,...) ##
More details to come
### fingerprint.py (Fingerprint, FingerprintMatrix, FingerprintSet) ###
### xios.py (XIOS, XiosEdge, MotifDB, Edge, Gspan) ###
### topology.py (Topology, SerialRNA, PairRNA, RNAstructure, Stem) ###




