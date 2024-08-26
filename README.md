# RNA XIOS Fingerprint

## Installation
1. Create a directory for the package
```
mkdir rna
cd rna
```
2. Install the XIOS/fingerprint package by cloning the git trepository. 
This will create a directory RNA with the XIOS/fingerprint code
```commandline
git clone https://github.com/RNA.git
```
2. Install Matthews RNAstructure package (command line version).
See https://rna.urmc.rochester.edu/RNAstructureDownload.html for more information. This will create a directory called RNAstructure. 
```
wget https://rna.urmc.rochester.edu/Releases/current/RNAstructureLinuxTextInterfaces64bit.tgz
tar -xvzf RNAstructureLinuxTextInterfaces64bit.tgz
```

## Making the motif library

## Constructing fingerprints from sequence
* Sequences should be in FASTA format and should contain only RNA bases ACGU
* Code locations below are with respect to the ***RNA*** directory.

Multiple steps can be run using *manager_new.py*, which manages submitting multiple jobs 
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

#### stochastic_to_xios.py
Calculate XIOS graph from sampled ct files ####
 * input: *.ct
 * output: *.xios
 * options
   - -c \<int> minimum number of counts to include paired bases as a stem (default=50)
   - -m \<int> minimum number of paired bases in a stem (default=3)




