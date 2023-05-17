#!/bin/bash
#SBATCH --job-name manager
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --account standby
#SBATCH --nodes=1
#SBATCH --ntasks=5

cd $SLURM_SUBMIT_DIR
python -V
python manager.py -r -w /scratch/bell/mgribsko/rna/RNA/manager_new/workflow0.yaml test0
