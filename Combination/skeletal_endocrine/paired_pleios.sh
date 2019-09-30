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

        sqlite3 /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/InfecPleiotropies_Apr2019/data/RAW_PROJECT/tests/PairPleiotropies/results/GWASpleiotropies.sqlite "SELECT DISTINCT SNPA,DiseaseA,RiskAllA,OnsetA,POS1,SNPB,DiseaseB,RiskAllB,OnsetB,POS2,ID,CHR FROM filteredPairs
WHERE R2 >= 0.8 AND ((GroupA IN ('skeletal system disease') AND GroupB IN ('endocrine system disease')) OR (GroupB IN ('skeletal system disease') AND GroupA IN ('endocrine system disease'))) 
AND CHR != '' ;" > /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Sep2019/results/Combination/skeletal_endocrine/pleios.txt

        ### Construct haplotypes with plink
        cat /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Sep2019/results/Combination/skeletal_endocrine/pleios.txt | while read line; do
         snpA=$(echo "$line" | cut -f 1);
         chr=$(echo "$line" | cut -f 12);
         echo $snpA;
         snpB=$(echo "$line" | cut -f 6);
         echo -e '*' ${snpA}'\t'${snpB} > /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Sep2019/results/Combination/skeletal_endocrine/snps.hlist;
         sed -i -e 's/ /\t/g' /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Sep2019/results/Combination/skeletal_endocrine/snps.hlist;
         /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/InfecPleiotropies_Apr2019/data/RAW_PROJECT/tests/PairPleiotropies/bin/plink-1.07-x86_64/plink --file /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/InfecPleiotropies_Apr2019/data/RAW_PROJECT/tests/PairPleiotropies/data/VCF_1000G/chr${chr}_CEU_genotypes --hap /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Sep2019/results/Combination/skeletal_endocrine/snps.hlist --hap-freq --noweb --out /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Sep2019/results/Combination/skeletal_endocrine/${snpA}_${snpB};
        grep -v LOCUS /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Sep2019/results/Combination/skeletal_endocrine/${snpA}_${snpB}.frq.hap | awk '{print $2,$3}' > /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Sep2019/results/Combination/skeletal_endocrine/${snpA}_${snpB}.fhtp
done

	### 3 #### Create pleiotropies agon/antagon with early/late classif
        cat /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Sep2019/results/Combination/skeletal_endocrine/pleios.txt | while read line; do
 snpA=$(echo "$line" | cut -f 1 );
        snpB=$(echo "$line" | cut -f 6 );
        nameA=$(echo "$line" | cut -f 2 );
        nameB=$(echo "$line" | cut -f 7 );
        RiskAllA=$(echo "$line" | cut -f 3 );
        RiskAllB=$(echo "$line" | cut -f 8 );
        PosA=$(echo "$line" | cut -f 5 );
        PosB=$(echo "$line" | cut -f 10 );
        OnsetA=$(echo "$line" | cut -f 4 );
        OnsetB=$(echo "$line" | cut -f 9 );
        chr=$(echo "$line" | cut -f 12 );
        riskHap=$(echo "$line" | awk -F"\t" '{print $3$8}');
### Python script to count pleiotropies in haplotypes
        for n in {10..60};
        do mkdir -p /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Sep2019/results/Combination/skeletal_endocrine/Age_threeshold_${n}/;
        python3 /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Sep2019/src/countPleiotropies.py /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Sep2019/results/Combination/skeletal_endocrine/${snpA}_${snpB}.fhtp ${n} ${riskHap} ${snpA} "${nameA}" ${RiskAllA} ${PosA} ${OnsetA} ${snpB} "${nameB}" ${RiskAllB} ${PosB} ${OnsetB} ${chr} /homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Sep2019/results/Combination/skeletal_endocrine/Age_threeshold_${n}/;
    done; done;
