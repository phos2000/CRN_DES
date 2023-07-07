#!/bin/bash

#SBATCH --mail-user=hanxuan.yu@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --mem=12G
#SBATCH --output=status/job%a.out
#SBATCH --array=1-100

module purge
module load GCC/10.2.0  OpenMPI/4.0.5 R-bundle-Bioconductor

Rscript --no-save sim_accre.R ${SLURM_ARRAY_TASK_ID}