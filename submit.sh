#!/bin/bash -l
#SBATCH --job-name=crc_nextflow
#SBATCH --output=nextflow_%j.log
#SBATCH --error=nextflow_%j.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=cpu

# Load Nextflow
module load nextflow/25.10.0-gcc-13.2.0

# Run pipeline
nextflow run main.nf -profile slurm -resume
