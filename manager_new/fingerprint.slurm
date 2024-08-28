#!/bin/bash
#SBATCH --job-name test_manager
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --account standby
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=4:00:00

module load anaconda
conda activate rna

cd $SLURM_SUBMIT_DIR
export set DATAPATH=/scratch/bell/mgribsko/rna/RNAstructure/data_tables
python ../RNA/manager_new/manager.py -w fingerprint.yaml -j 32 testmanager


