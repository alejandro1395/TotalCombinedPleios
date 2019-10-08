#!/bin/bash

#set the job name
#SBATCH --job-name=Plots_loop
#SBATCH --cpus-per-per-task = 1
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --time=01:00:00

#run the application
#VARIABLES
##DISEASE GROUPS

skin="skin disease"
skeletal="skeletal system disease"
declare -a respiratory=("respiratory system disease" "Abnormality of the respiratory system")
reproductive="reproductive system disease"
declare -a cancer=("neoplasm" "Meningioma")
immune="immune system disease"
infectious="infectious disease"
declare -a kidney=("kidney disease" "rapid kidney function decline")
declare -a cardio=("cardiovascular disease" "Abnormality of the cardiovascular system")
declare -a metabolic=("digestive system disease" "metabolic disease" "Abnormality of metabolism/homeostasis")
declare -a headneck=("Abnormality of head or neck" "head and neck disorder")
declare -a nervous=("Abnormality of the nervous system" "nervous system disease")
eye="eye disease"
endocrine="endocrine system disease"


#output folder
OUTPUT=/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Sep2019/results/Combination/
touch ${OUTPUT}freqs_Agon_Antagon.tsv
echo -e "Combined_diseases\tAgonistic\tAntagonistic" >> ${OUTPUT}freqs_Agon_Antagon.tsv


#LOOP THROUGH ALL GROUPS OF DISEASES ALONE

groups=("skin" "skeletal" "respiratory" "reproductive" "infectious" "immune" "kidney" "cardio" "metabolic" "headneck" "nervous" "eye" "endocrine" "cancer")
max=${#groups[@]}
for ((idxA=0; idxA<max; idxA++)); do
for ((idxB=idxA; idxB<max; idxB++)); do
         # iterate idxB from idxA to length

        j="${groups[$idxA]}"
        i="$j[@]"
        a=("${!i}")
        g="${groups[$idxB]}"
        p="$g[@]"
        b=("${!p}")
        echo ${j}_${g}

	INPUT=/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Sep2019/results/Combination/${j}_${g}/
#Loop over all ages Single Pleios

Early_Late_Agonistic_Single=$(cat ${INPUT}Age_threeshold_10/SinglePleiosLoc.tsv | awk -F $'\t' 'BEGIN {OFS = FS} ($8=="Early-Late Agonistic"){++count} END {print count}')
Late_Early_Agonistic_Single=$(cat ${INPUT}Age_threeshold_10/SinglePleiosLoc.tsv | awk -F $'\t' 'BEGIN {OFS = FS} ($8=="Late-Early Agonistic"){++count} END {print count}')
Late_Late_Agonistic_Single=$(cat ${INPUT}Age_threeshold_10/SinglePleiosLoc.tsv | awk -F $'\t' 'BEGIN {OFS = FS} ($8=="Late-Late Agonistic"){++count} END {print count}')
Early_Early_Agonistic_Single=$(cat ${INPUT}Age_threeshold_10/SinglePleiosLoc.tsv | awk -F $'\t' 'BEGIN {OFS = FS} ($8=="Early-Early Agonistic"){++count} END {print count}')

if [ -z "${Early_Late_Agonistic_Single}" ]; then
    Early_Late_Agonistic_Single=0
fi
if [ -z "${Late_Early_Agonistic_Single}" ]; then
    Late_Early_Agonistic_Single=0
fi
if [ -z "${Late_Late_Agonistic_Single}" ]; then
    Late_Late_Agonistic_Single=0
fi
if [ -z "${Early_Early_Agonistic_Single}" ]; then
    Early_Early_Agonistic_Single=0
fi


#antagonistic singles

Early_Late_Antagonistic_Single=$(cat ${INPUT}Age_threeshold_10/SinglePleiosLoc.tsv | awk -F $'\t' 'BEGIN {OFS = FS} ($8=="Early-Late Antagonistic"){++count} END {print count}')
Late_Early_Antagonistic_Single=$(cat ${INPUT}Age_threeshold_10/SinglePleiosLoc.tsv | awk -F $'\t' 'BEGIN {OFS = FS} ($8=="Late-Early Antagonistic"){++count} END {print count}')
Late_Late_Antagonistic_Single=$(cat ${INPUT}Age_threeshold_10/SinglePleiosLoc.tsv | awk -F $'\t' 'BEGIN {OFS = FS} ($8=="Late-Late Antagonistic"){++count} END {print count}')
Early_Early_Antagonistic_Single=$(cat ${INPUT}Age_threeshold_10/SinglePleiosLoc.tsv | awk -F $'\t' 'BEGIN {OFS = FS} ($8=="Early-Early Antagonistic"){++count} END {print count}')

if [ -z "${Early_Late_Antagonistic_Single}" ]; then
    Early_Late_Antagonistic_Single=0
fi
if [ -z "${Late_Early_Antagonistic_Single}" ]; then
    Late_Early_Antagonistic_Single=0
fi
if [ -z "${Late_Late_Antagonistic_Single}" ]; then
    Late_Late_Antagonistic_Single=0
fi
if [ -z "${Early_Early_Antagonistic_Single}" ]; then
    Early_Early_Antagonistic_Single=0
fi

#agonistic pairs

Early_Late_Agonistic_Pair=$(wc -l ${INPUT}Age_threeshold_10/Agon_early_late | cut -f 1 -d " ")
Early_Early_Agonistic_Pair=$(wc -l ${INPUT}Age_threeshold_10/Agon_early_early | cut -f 1 -d " ")
Late_Late_Agonistic_Pair=$(wc -l ${INPUT}Age_threeshold_10/Agon_late_late | cut -f 1 -d " ")


if [ -z "${Early_Late_Agonistic_Pair}" ]; then
    Early_Late_Agonistic_Pair=0
fi
if [ -z "${Early_Early_Agonistic_Pair}" ]; then
    Early_Early_Agonistic_Pair=0
fi
if [ -z "${Late_Late_Agonistic_Pair}" ]; then
    Late_Late_Agonistic_Pair=0
fi

#Antagonistic pairs

Early_Late_Antagonistic_Pair=$(wc -l ${INPUT}Age_threeshold_10/Antagon_early_late | cut -f 1 -d " ")
Early_Early_Antagonistic_Pair=$(wc -l ${INPUT}Age_threeshold_10/Antagon_early_early | cut -f 1 -d " ")
Late_Late_Antagonistic_Pair=$(wc -l ${INPUT}Age_threeshold_10/Antagon_late_late | cut -f 1 -d " ")


if [ -z "${Early_Late_Antagonistic_Pair}" ]; then
    Early_Late_Antagonistic_Pair=0
fi
if [ -z "${Early_Early_Antagonistic_Pair}" ]; then
    Early_Early_Antagonistic_Pair=0
fi
if [ -z "${Late_Late_Antagonistic_Pair}" ]; then
    Late_Late_Antagonistic_Pair=0
fi


#Print sums into output
Agonistic=$((${Early_Late_Agonistic_Pair} + ${Early_Late_Agonistic_Single} + ${Late_Early_Agonistic_Single} \
+ ${Late_Late_Agonistic_Pair} + ${Late_Late_Agonistic_Single} + ${Early_Early_Agonistic_Single} + ${Early_Early_Agonistic_Pair}))
Antagonistic=$((${Early_Late_Antagonistic_Pair} + ${Early_Late_Antagonistic_Single} + ${Late_Early_Antagonistic_Single} \
+ ${Late_Late_Antagonistic_Pair} + ${Late_Late_Antagonistic_Single} + ${Early_Early_Antagonistic_Single} + ${Early_Early_Antagonistic_Pair}))

echo -e "${j}_${g}\t${Agonistic}\t${Antagonistic}" >> ${OUTPUT}freqs_Agon_Antagon.tsv
done;done;
