#!/usr/bin/python3

import pandas as pd
import numpy as np
import requests, sys
from itertools import combinations
#import seaborn as sns
from scipy import stats
import matplotlib.pyplot as plt
import pickle
from collections import Counter
from matplotlib.backends.backend_pdf import PdfPages
import copy
from scipy.stats import sem, t
from scipy import mean
import re
import os

"""
CREATING DATABASE WITHDIFF SINGLE PLEIOTROPIES
"""

#VARIABLES
PATH = "/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Jul2019/Cancer/"
OUTPATH = sys.argv[2]
catalog = pd.read_csv(PATH + "data/GWAS_Age_merged.csv", sep='\t', low_memory=False)#panda creation
LD = pd.read_csv("~/scratch/PhD_EvoGenomics/1st_year/InfecPleiotropies_Apr2019/data/RAW_PROJECT/CEU_mod.csv", sep='\t')
SNPs_multiple_diseases = {}
out_file1 = OUTPATH + "SinglePleiotropiesOut.tsv"
out_file2 = OUTPATH + "SinglePleiosLoc.tsv"
age_threeshold = int(sys.argv[1])
diseaseA_groups = eval(sys.argv[3])
if type(diseaseA_groups) is str:
    diseaseA_groups = [str(diseaseA_groups)]
diseaseB_groups = eval(sys.argv[4])
if type(diseaseB_groups) is str:
    diseaseB_groups = [str(diseaseB_groups)]

## 1st strategy -> Single SNPs influencing more than one disease

def create_diseases_dictionary_for_each_SNP(catalog_fun, SNPs_dict):
    for i, row in catalog_fun.iterrows():
        if row['SNPS'] not in SNPs_dict.keys():
            SNPs_dict[row['SNPS']] = [row['DISEASE/TRAIT']]
        else:
            if row['DISEASE/TRAIT'] not in SNPs_dict[row['SNPS']]:
                SNPs_dict[row['SNPS']].append(row['DISEASE/TRAIT'])
    return SNPs_dict


def filter_pleiotropic_single_SNPs(SNPs_dict):
    pleiotropic_singles = {k:v for k,v in SNPs_dict.items() if len(SNPs_dict[k])>=2}
    return pleiotropic_singles

def get_subset_disease_Pleiotropies(catalog_fun, pleio_dict, diseasesA, diseasesB):
    disease = []
    pleios = pd.DataFrame(columns=['SNP', 'Diseases', 'Group',
    'REGION', 'CHR_ID', 'CHR_POS', 'PleioType'])
    for entry in pleio_dict:
        subset = catalog_fun.loc[catalog_fun['SNPS'] == entry]
        subset["PleioType"] = ""
        if any(str(elem) in diseasesA for elem in subset['Group'].tolist()) and any(str(elem) in diseasesB for elem in subset['Group'].tolist()):
            diseasesA_subset = subset[subset['Group'].isin(diseasesA)]
            diseasesB_subset = subset[subset['Group'].isin(diseasesB)]
            diseasesA_riskalleles = retrieve_subset_of_column_values(diseasesA_subset, "Risk_allele")
            diseasesB_riskalleles = retrieve_subset_of_column_values(diseasesB_subset, "Risk_allele")
            diseasesA_ages = retrieve_subset_of_column_values(diseasesA_subset, "Onset")
            diseasesB_ages = retrieve_subset_of_column_values(diseasesB_subset, "Onset")
            if set(diseasesA_riskalleles) == set(diseasesB_riskalleles):
                type = "Agonistic"
            else:
                type = "Antagonistic"
            pleio_type = define_early_late(age_threeshold, diseasesB_ages, diseasesB_ages, type)
            pleios = pleios.append({'SNP': entry, 'Diseases': set(subset['Disease'].tolist()),
            'Group': set(subset['Group'].tolist()), 'REGION': set(subset['Disease'].tolist()),
            'CHR_ID': set(subset['CHR_ID'].tolist()), 'CHR_POS': set(subset['CHR_POS'].tolist()),
            'PleioType': pleio_type}, ignore_index=True)
    return pleios


def retrieve_subset_of_column_values(dataset, column_name):
    dataset_values = dataset[column_name].tolist()
    return dataset_values


def define_early_late(threeshold, dis1age, dis2age, type):
    if all(x <= threeshold for x in dis1age) and all(x <= threeshold for x in dis2age):
        pleiotype = "Early-Early " + type
    elif all(x > threeshold for x in dis1age) and all(x > threeshold for x in dis2age):
        pleiotype = "Late-Late " + type
    elif any(x <= threeshold for x in dis1age) and any(x > threeshold for x in dis2age):
        pleiotype = "Early-Late " + type
    elif any(x > threeshold for x in dis1age) and any(x <= threeshold for x in dis2age):
        pleiotype = "Late-Early " + type
    return pleiotype


"""
def get_subset_Cancer_Pleiotropies(table, array):
    for i, row in table.iterrows():
        if row['Group'] == "neoplasm" or row['Group'] == "Meningioma":
            array.append(row['SNP'])
    table_subset1 = table[table.SNP.isin(array)]
    table_subset1 = table_subset1.sort_values(by=['SNP'])
    return table_subset1

def get_subset_immunitary_Pleiotropies(table, array):
    for i, row in table.iterrows():
        if row['Group'] == "immune system disease":
            array.append(row['SNP'])
    table_subset2 = table[table.SNP.isin(array)]
    table_subset2 = table_subset2.sort_values(by=['SNP'])
    return table_subset2

def agonistic_or_antagonistic(table):
    Agonistic = 0
    Antagonistic = 0
    for snp in table.SNP.unique().tolist():
        df = table.loc[table['SNP'] == snp]
        if len(df.RiskAll.unique().tolist()) == 1:
            Agonistic += 1
        else:
            Antagonistic +=1
    #while table.loc[df['column_name'] == some_value]
    return Agonistic, Antagonistic
"""


#def sadasfsafsaf():

#MAIN

SNPs_multiple_diseases = create_diseases_dictionary_for_each_SNP(catalog, SNPs_multiple_diseases)
#print(SNPs_multiple_diseases)
pleiotropic_singles = filter_pleiotropic_single_SNPs(SNPs_multiple_diseases)
#for key in pleiotropic_singles.keys():
    #print(key, pleiotropic_singles[key])
pleios_diseasesAB = get_subset_disease_Pleiotropies(catalog, pleiotropic_singles, diseaseA_groups, diseaseB_groups)
pleios_diseasesAB.to_csv(out_file2, sep='\t')

"""
SinglePleio_subset1 = get_subset_infectious_Pleiotropies(SinglePleio, SingleInfectSNPs)
SinglePleio_subset2 = get_subset_immunitary_Pleiotropies(SinglePleio, SingleInfectSNPs)
SinglePleio_subset2.to_csv(out_file, sep='\t')

Agonistic, Antagonistic = agonistic_or_antagonistic(SinglePleio_subset2)
print(Agonistic, Antagonistic)
"""
