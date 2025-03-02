#!/bin/bash
#SBATCH --job-name p_310_15_5
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --account=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=4:00:00

cd $SLURM_SUBMIT_DIR
export set DATAPATH=/scratch/scholar/ashome/rna/RNAstructure/data_tables/
python /scratch/scholar/ashome/rna/RNA/manager_new/manager.py -w /scratch/scholar/ashome/rna/RNA/manager_new/workflow_new.yaml -j 32 p_310_15_5
