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
* Sequences should be in FASTA format and should contain only 