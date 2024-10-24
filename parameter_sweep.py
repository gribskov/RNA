"""=================================================================================================
create workflow files and slurm files for testing parameters of xios file generation pipeline using
RNAstructure:Partition
RNAstructure:Stochastic
stochastic_to_xios.py

parameters are set below in the main program
T = folding temperature: 
F = minimum occurrence of base-pairs in stochastic sample:
S = minimum stem length:

output is produced in the current directory when the program is run. 
the base variable must be defined in main to point to the directory where RNA and RNAstructure
are installed. currently base = '/scratch/bell/mgribsko/rna/'

All parameter settings use the same random seed in stochastic. This seed depends on the current
time so it will be different for every run of parameter_sweep

------------------
Running the jobs
All jobs can be run using a simple bash script such as

for f in *.slurm; do
  echo "$f"
  sbatch $f
done
================================================================================================="""
from itertools import product
from time import time, sleep


def slurm_template():
"""-------------------------------------------------------------------------------------------------
return a string with the header for the Slurm job file. XPROJECT, XACCOUNT, XXIOS_EXE are replaced
with real values in main
-------------------------------------------------------------------------------------------------"""
    return """#!/bin/bash
#SBATCH --job-name XPROJECT
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --account XACCOUNT
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=4:00:00

cd $SLURM_SUBMIT_DIR
export set DATAPATH=XRNASTRUCTURE_DATA
python XXIOS_EXEmanager_new/manager.py -w XWORKFLOW -j 32 XPROJECT
"""


def workflow_template():
"""-------------------------------------------------------------------------------------------------
return a string with the yaml formatted workflow. XRNASTRUCTURE, XXIOS_EXE, and XPROJECT are
replaced with real values in main
-------------------------------------------------------------------------------------------------"""
    return """
---
  definitions:
#   code locations
    python: python
    base: /scratch/bell/mgribsko/rna
    RNAstructure: XRNASTRUCTURE_EXE
    xiosexe: XXIOS_EXE
#   project input/output directories
    fasta: $base/RNA/data/curated_fasta
    project: XPROJECT
    partition: $project/partition
    stochastic: $project/stochastic
    xios: $project/xios
  stage:
    partition:
      command: $RNAstructure/partition  $in $out $option
      option:
        -t XT
        -q
      in: $fasta/*.fa
      out: :$partition/%in.replace('.fa', '.pfs')
    stochastic:
      command: $RNAstructure/stochastic $in $out $option
      option: -s XR
      in: $partition/*.pfs
      out: :$stochastic/%in.replace('.pfs', '.ct')
    xios:
      command: $python $xiosexe/stochastic_to_xios.py $option $in $out
      option: -m XS -c XF
      in: $stochastic/*.ct
      out: :$xios/%in.replace('.ct', '.xios')
"""


"""=================================================================================================
# main program

base = root directory contain RNA and RNAStructure code (/scratch/bell/mgribsko/rna/)
XACCOUNT: Slurm account to run under (standby)
XRNASTRUCTURE_DATA': RNAstructure energy tables directory (f'{base}RNAstructure/data_tables/')
XRNASTRUCTURE_EXE': RNAstructure program executables (f'{base}RNAstructure/exe/')
XXIOS_EXE': XIOS program executables (f'{base}RNA/'})

subs['XT'] = list of folding temperatures to survey (partition)
subs['XF'] = list of minimum fraction to survey (stochastic)
subs['XR'] = random number seed for sampling from ensemble (stochastic), currently int(time()/101)
subs['XS'] = list of minimum stem length to survey (stochastic_to_xios.py)
subs['XPROJECT'] = name for parameter setting directory, currently f'p_{t}_{f}_{s}'
subs['XWORKFLOW'] = name for workflow file, currently f'{subs["XPROJECT"]}.workflow'
================================================================================================="""
if __name__ == '__main__':
    T = [260 + i * 10 for i in range(9)]
    F = [10, 20, 30, 40, 50, 60]
    S = [2, 3, 4]

    base = '/scratch/bell/mgribsko/rna/'
    subs = {'XACCOUNT': 'standby',
            'XRNASTRUCTURE_DATA': f'{base}RNAstructure/data_tables/',
            'XRNASTRUCTURE_EXE': f'{base}RNAstructure/exe/',
            'XXIOS_EXE': f'{base}RNA/'}

    for t, f, s in product(T, F, S):
        print(f'{t}  {f}  {s}')
        subs['XT'] = f'{t}'
        subs['XF'] = f'{f}'
        subs['XS'] = f'{s}'
        r = int(time()/101)
        if not r % 2:
            r += 199
        subs['XR'] = f'{r}'

        subs['XPROJECT'] = f'p_{t}_{f}_{s}'
        subs['XWORKFLOW'] = f'{subs["XPROJECT"]}.workflow'
        slurm = slurm_template()
        workflow = workflow_template()
        for substitution in subs:
            slurm = slurm.replace(substitution, subs[substitution])
            workflow = workflow.replace(substitution, subs[substitution])

        slurmout = open(f'{subs["XPROJECT"]}.slurm', 'w')
        slurmout.write(slurm)
        slurmout.close()        
        print(slurm)

        workout = open(subs['XWORKFLOW'], 'w')
        workout.write(workflow)
        workout.close()
        print(workflow)
        sleep(0.01)
