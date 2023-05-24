"""=================================================================================================
create workflow files and slurm files for testing parameters of xios file generation pipeline using
RNAstructure:Partition
RNAstructure:Stochastic
stochastic_to_xios.py

parameters are
T = folding temperature: 280 - 340 by 20
F = minimum occurrence of base-pairs in stochastic sample: 30, 50, 80
S = minimum stem length: 3, 4
================================================================================================="""
from itertools import product
from time import time, sleep


def slurm_template():
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
    return """
---
  definitions:
#   code locations
    python: python
    RNAstructure: XRNASTRUCTURE_EXE
    xiosexe: XXIOS_EXE
#   project input/output directories
    fasta: ../standards/fasta_fixed
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
================================================================================================="""
if __name__ == '__main__':
    T = [280 + i * 10 for i in range(7)]
    F = [30, 50, 80]
    S = [3, 4]

    subs = {'XACCOUNT': 'standby',
            'XRNASTRUCTURE_DATA': '~/RNAstructure/data_tables/',
            'XRNASTRUCTURE_EXE': '~/RNAstructure/exe/',
            'XXIOS_EXE': '~/PycharmProjects/RNA/'}

    for t, f, s in product(T, F, S):
        print(f'{t}  {f}  {s}')
        subs['XT'] = f'{t}'
        subs['XF'] = f'{f}'
        subs['XS'] = f'{s}'
        r = int(time()/1001)
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
