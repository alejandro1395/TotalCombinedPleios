#!/bin/bash

#set the job name
#SBATCH --job-name=Haps_Onset_Immune
#SBATCH --cpus-per-per-task = 1
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --time=01:00:00

module load Python

### Python script to count pleiotropies in haplotypes
	for t in {10..60};
  	 do mkdir -p /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Sep2019/results/Alone/respiratory/Age_threeshold_${t}/;	
\python3 /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Sep2019/src/SinglePleios_allGroupsAlone.py ${t} /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Sep2019/results/Alone/respiratory/Age_threeshold_${t}/ "('respiratory system disease', 'Abnormality of the respiratory system')" 
    done;
