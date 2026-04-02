#!/bin/bash

#SBATCH --array=1 
#SBATCH --mem=50g
#SBATCH --time=140:00:00

module load Conda/3
conda activate py_env

python fgDNN_ADHD.py $SLURM_ARRAY_TASK_ID