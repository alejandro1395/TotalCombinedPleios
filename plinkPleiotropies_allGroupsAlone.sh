#!/bin/bash

#set the job name
#SBATCH --job-name=Haps_Onset_Cancer
#SBATCH --cpus-per-per-task = 1
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --time=03:00:00

module load Python

#RUN HE APPLICATION
#PATHS
OUTPUT=/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Sep2019/results/
SRC=/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Sep2019/src/
INPUT=/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/InfecPleiotropies_Apr2019/data/RAW_PROJECT/tests/PairPleiotropies/results/
VCF=/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/InfecPleiotropies_Apr2019/data/RAW_PROJECT/tests/PairPleiotropies/data/VCF_1000G/
BIN=/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/InfecPleiotropies_Apr2019/data/RAW_PROJECT/tests/PairPleiotropies/bin/plink-1.07-x86_64/
mkdir -p ${OUTPUT}Alone
mkdir -p ${OUTPUT}Combination
mkdir -p ${SRC}Alone
mkdir -p ${SRC}Combination


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

#LOOP THROUGH ALL GROUPS OF DISEASES ALONE

groups=("skin" "skeletal" "respiratory" "reproductive" "infectious" "immune" "kidney" "cardio" "metabolic" "headneck" "nervous" "eye" "endocrine" "cancer")
for i in "${groups[@]}"; do
    j="$i[@]"
    a=("${!j}")
	mkdir -p ${SRC}Alone/${i}/
	

#FORMATING GROUPS OF DISEASES IN STRINGS
	separator="', '"

	diseasesGroup="$( printf "${separator}%s" "${a[@]}" )"
	mkdir -p ${OUTPUT}Alone/${i}
	diseasesGroup="${diseasesGroup:${#separator}}" # remove leading separator
	diseases=$(echo "('$diseasesGroup')")
	
	FILE=${SRC}Alone/${i}/paired_pleios.sh
	
cat  > $FILE << EOT 
#!/bin/bash
#set the job name
#SBATCH --job-name=Haps_Onset_Cancer
#SBATCH --cpus-per-per-task = 1
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --time=03:00:00

module load Python
	
#### Crear haplotipos con Plink. Solo lo hace una vez por valor de R2
### Esta comanda de aqui abajo selecciona pleiotropias por el valor de r2 que tienen en la tabla.
	
	sqlite3 ${INPUT}GWASpleiotropies.sqlite "SELECT DISTINCT SNPA,DiseaseA,RiskAllA,OnsetA,POS1,SNPB,DiseaseB,RiskAllB,OnsetB,POS2,ID,CHR FROM filteredPairs 
WHERE R2 >= 0.8 AND (GroupA IN $(echo $diseases) OR GroupB IN $(echo $diseases)) AND CHR != '' ;" > ${OUTPUT}Alone/${i}/pleios.txt

	### Construct haplotypes with plink
	cat ${OUTPUT}Alone/${i}/pleios.txt | while read line; do
   	 snpA=\$(echo "\$line" | cut -f 1);
   	 chr=\$(echo "\$line" | cut -f 12);
   	 echo \$snpA;
   	 snpB=\$(echo "\$line" | cut -f 6);
   	 echo -e '*' \${snpA}'\t'\${snpB} > ${OUTPUT}Alone/${i}/snps.hlist;
   	 sed -i -e 's/ /\t/g' ${OUTPUT}Alone/${i}/snps.hlist;
   	 ${BIN}plink --file ${VCF}chr\${chr}_CEU_genotypes \
--hap ${OUTPUT}Alone/${i}/snps.hlist \
--hap-freq \
--noweb \
--out ${OUTPUT}Alone/${i}/\${snpA}_\${snpB};
	grep -v LOCUS ${OUTPUT}Alone/${i}/\${snpA}_\${snpB}.frq.hap | awk '{print \$2,\$3}' > ${OUTPUT}Alone/${i}/\${snpA}_\${snpB}.fhtp
	done

	### 3 #### Create pleiotropies agon/antagon with early/late classif
	cat ${OUTPUT}Alone/${i}/pleios.txt | while read line; do
 snpA=\$(echo "\$line" | cut -f 1 );
        snpB=\$(echo "\$line" | cut -f 6 );
        nameA=\$(echo "\$line" | cut -f 2 );
        nameB=\$(echo "\$line" | cut -f 7 );
        RiskAllA=\$(echo "\$line" | cut -f 3 );
        RiskAllB=\$(echo "\$line" | cut -f 8 );
        PosA=\$(echo "\$line" | cut -f 5 );
        PosB=\$(echo "\$line" | cut -f 10 );
        OnsetA=\$(echo "\$line" | cut -f 4 );
        OnsetB=\$(echo "\$line" | cut -f 9 );
        chr=\$(echo "\$line" | cut -f 12 );
        riskHap=\$(echo "\$line" | awk -F"\t" '{print \$3\$8}');
### Python script to count pleiotropies in haplotypes
        for i in {10..60};
        do mkdir -p ${OUTPUT}Alone/${i}/Age_threeshold_\${i}/;
        python3 ./countPleiotropies.py ${OUTPUT}Alone/${i}/\${snpA}_\${snpB}.fhtp \${i} \${riskHap} \${snpA} "\${nameA}" \${RiskAllA} \${PosA} \
\${OnsetA} \${snpB} "\${nameB}" \${RiskAllB} \${PosB} \${OnsetB} \${chr} ${OUTPUT}Alone/${i}/Age_threeshold_\${i}/;
    done; done;
EOT
	sbatch ${SRC}Alone/${i}/paired_pleios.sh
cat > ${SRC}Alone/${i}/SinglePleios.sh << EOT	
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
  	 do mkdir -p ${OUTPUT}Alone/${i}/Age_threeshold_\${t}/;	
\python3 ${SRC}SinglePleios_allGroupsAlone.py \${t} ${OUTPUT}Alone/${i}/Age_threeshold_\${t}/ "$diseases" 
    done;
EOT
	sbatch ${SRC}Alone/${i}/SinglePleios.sh 
done
