#!/bin/bash
#SBATCH --account=egitim
#SBATCH --partition=barbun
#SBATCH --ntasks-per-node=4
#SBATCH --output=logs/slurm/FASTQC_BEFORE_TRIMMING.out
#SBATCH --error=logs/slurm/FASTQC_BEFORE_TRIMMING.err

PATH=${PATH}:/truba/home/egitim/miniconda3/envs/aMeta/bin/

mkdir -p results/FASTQC_BEFORE_TRIMMING logs/FASTQC_BEFORE_TRIMMING
for sample in $(ls data/*.fastq.gz); do
  sample_name=$(basename $sample .fastq.gz)
  fastqc $sample --threads 4 --nogroup --outdir results/FASTQC_BEFORE_TRIMMING &> logs/FASTQC_BEFORE_TRIMMING/$sample_name.log;
done