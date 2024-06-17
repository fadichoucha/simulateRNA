#!/bin/bash
#SBATCH --job-name=rnaseq_simulation
#SBATCH --output=../output/rnaseq_simulation.out
#SBATCH --error=../output/rnaseq_simulation.err
#SBATCH --time=24:00:00
#SBATCH --partition=your_partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G


module load anaconda/3
source activate rnaseq_simulation

Rscript rnaseq_simulation.R

conda deactivate
