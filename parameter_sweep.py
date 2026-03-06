"""=================================================================================================
create workflow files and slurm files for testing parameters of XIOS/fingerprint  generation pipeline using
    RNAstructure:Partition
    RNAstructure:Stochastic
    stochastic_to_xios.py
    fingerprint_random.py
data/parameter_sweep.workflow provides a template command that this script updates with desired parameters

parameters are set below in the main program
T = folding temperature: 
F = minimum occurrence of base-pairs in stochastic sample:
S = minimum stem length:

output is determined by the project in the workflow, which is set by this script
paths to data and executable are also defined in workflow

RNAstructure:stochastic requires a random seed. this is set based on current time so repeat runs will no necessarily
produce exactly the same result
TODO mechanism to exactly replicate by providing a defined seed
================================================================================================="""
from itertools import product
from time import time, sleep
import yaml
from manager_new.manager import Command

"""-------------------------------------------------------------------------------------------------
Template for the Slurm job file header. Most of the information is filled in main:
XPROJECT
XACCOUNT
XCPUS
XENVIRONMENT
XRNASTRUCTUREDATA
XMANAGER
-------------------------------------------------------------------------------------------------"""
slurm_template = """
    #!/bin/bash
    #SBATCH --job-name <XPROJECT>
    #SBATCH --output=%x_%j.out
    #SBATCH --error=%x_%j.err
    #SBATCH --account <XACCOUNT>
    #SBATCH --nodes=1
    #SBATCH --ntasks=<XCPUS>
    #SBATCH --time=4:00:00
    
    module load anaconda
    conda activate <XENVIRONMENT>
    
    export set DATAPATH=<XRNASTRUCTURE_DATA>
    python <XMANAGER> <XPROJECT>.workflow -j <XCPUS>
"""

"""=================================================================================================
# main program

slurm_template variables provided by main:
XPROJECT            project name for this parameter set; defines name for workflow file    
XACCOUNT            account for slurm job
XCPUS               number of CPUs to request through slurm; number of process to keep running (manager.py)
XENVIRONMENT        conda environment for python (stochastic_to_xios.py and fingerprint_random.py)
XRNASTRUCTUREDATA   path to RNAstructure energy tables
XMANAGER            path to XIOS workflow manager

variables substituted into workflow template
subs['XT'] = list of folding temperatures to survey (partition)
subs['XF'] = list of minimum fraction to survey (stochastic)
subs['XR'] = random number seed for sampling from ensemble (stochastic), currently int(time()/101)
subs['XS'] = list of minimum stem length to survey (stochastic_to_xios.py)
subs['XPROJECT'] = name for output directory for this parameter set, currently f'p_{t}_{f}_{s}'
subs['XWORKFLOW'] = name for workflow file, currently f'{subs["XPROJECT"]}.workflow'
================================================================================================="""
if __name__ == '__main__':
    # Temperature, Fraction, Stem size
    T = [t for t in range(260, 360, 10)]
    F = [f for f in range(10, 225, 10)]
    S = [s for s in range(2, 6)]

    njobs = len(T) * len(F) * len(S)
    print(f'[arameter_sweep')
    print(f'\tTemperature: {T}')
    print(f'\tFraction: {F}')
    print(f'\tStem minimum size: {S}')
    print(f'\t Total jobs: {njobs}')

    base = '/scratch/bell/mgribsko/rna/'
    subs = {'XACCOUNT': 'standby',
            'XCPUS': '128',
            'XENVIRONMENT': 'RNA26',
            'XPROJECT':None,
            'XWORKFLOW': None,
            'XRNASTRUCTURE_DATA': f'{base}RNAstructure/data_tables/',
            'XRNASTRUCTURE_EXE': f'{base}RNAstructure/exe/',
            'XMANAGER': f'{base}RNA/manager_new/manager.py'}
    workflow_template = f'{base}/RNA/data/parameter_sweep.workflow'

    for t, f, s in product(T, F, S):
        print(f'{t}  {f}  {s}')
        subs['XT'] = f'{t}'
        subs['XF'] = f'{f}'
        subs['XM'] = f'{s}'


        subs['XPROJECT'] = f'p_{t}_{f}_{s}'
        subs['XWORKFLOW'] = f'{subs["XPROJECT"]}.workflow'
        r = int(time() / 101)
        if not r % 2:
            # avoid factor of 2 in random number seed
            r += 199
        subs['XSEED'] = f'{r}'

        # fill in placeholders in workflow template and save
        template = Command('')
        template = open('data/parameter_sweep.workflow','r').read()

        for substitution in subs:
            target = f'<{substitution}>'
            template = template.replace(target, subs[substitution])

        wf = open(f'{subs["XPROJECT"]}.workflow', 'w')
        wf.write(template)
        wf.close()

        # fill in placeholders in slurm job templates
        slurm = slurm_template
        for substitution in subs:
            target = f'<{substitution}>'
            slurm = slurm.replace(target, subs[substitution])

        slurmout = open(f'{subs["XPROJECT"]}.slurm', 'w')
        slurmout.write(slurm)
        slurmout.close()
        print(slurm)

        sleep(0.01)
