#!/bin/bash

#set the job name
#SBATCH --job-name=Haps_Onset_Single
#SBATCH --cpus-per-per-task = 1
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --time=01:00:00

module load Python

### Python script to count pleiotropies in haplotypes
        for m in {10..60};
         do mkdir -p /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Sep2019/results/Combination/metabolic_respiratory/Age_threeshold_${m}/;
\python3 /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Sep2019/src/SinglePleios_allGroupsCombinations.py ${m} /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Sep2019/results/Combination/metabolic_respiratory/Age_threeshold_${m}/ "('digestive system disease', 'metabolic disease', 'Abnormality of metabolism/homeostasis')" "('respiratory system disease', 'Abnormality of the respiratory system')"
    done;
