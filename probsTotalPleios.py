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
from time import sleep

#VARIABLES
Pairs_of_diseases = [["Cancer", "Immune"], ["Cancer", "Metabolic"], ["Metabolic", "Immune"]]
Pleiotropies=["Agon_early_early", "Agon_early_late", "Agon_early_early",
"Antagon_early_early", "Antagon_early_late", "Agon_early_early"]
OUTPATH="/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Jul2019/Combinations/plots/"
total_pleios = 1803
probs_pleios = []
total_comb_pleios = []

"""
Functions
"""

def count_number_of_pleios(PATH):
	if os.path.exists(PATH):
		with open(PATH, "r") as f:
			size=sum(1 for _ in f)
	else:
		size=0
	return(size)


def count_AP_EarlyLate_single(PATH):
	Singlesdb = pd.read_csv(PATH, sep="\t")
	len1 = 0
	len2 = 0
	len1 += len(Singlesdb[Singlesdb['PleioType'] == 'Early-Late Antagonistic'])
	len2 += len(Singlesdb[Singlesdb['PleioType'] == 'Late-Early Antagonistic'])
	len_tot = len1 + len2
	return(len_tot)



#LOOPS MAIN
#First of all, get frequency and multiply to know probability of AP in that gene
for pair in Pairs_of_diseases:
	pair_APs = []
	totals_disease = []
	for dis in pair:
		#Paired SNP pleiotropies
		INPUT_PATH="/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Jul2019/" + dis + "/results/Onset/"
		total_pair_disease = count_number_of_pleios(INPUT_PATH + dis + "_pleios/pleios.txt")
		#Single SNP pleiotropies
		total_single_disease = count_number_of_pleios(INPUT_PATH +
		"SinglePleios/Age_threeshold_10/SinglePleiosLoc.tsv") - 1
		total_disease = total_pair_disease + total_single_disease
		totals_disease.append(total_disease)
	#TOTAL PROB
	freq_disease1 = (totals_disease[0])/total_pleios
	freq_disease2 = (totals_disease[1])/total_pleios
	Prob_together = freq_disease1 * freq_disease2
	probs_pleios.append(Prob_together)



	#Now we go with the count of OBSERVED combined
	#First with the paired
	INPUT_COMB_PATH="/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Jul2019/Combinations/results/" + \
	pair[0] + pair[1] + "_pleios/pleios.txt"
	Comb_pair_val = count_number_of_pleios(INPUT_COMB_PATH)
	INPUT_Single_COMB_PATH = "/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Jul2019/Combinations/results/" + \
	"SinglePleios_" + pair[0] + pair[1] + "/Age_threeshold_10/SinglePleiosLoc.tsv"
	single_comb = count_number_of_pleios(INPUT_Single_COMB_PATH)
	Total_comb = Comb_pair_val + single_comb
	total_comb_pleios.append(Total_comb)

Prob_none = 1 - sum(probs_pleios)
Rest_comb = total_pleios - sum(total_comb_pleios)

#Finally, we print the ouput in a new folders
with open(OUTPATH + "freqs_total.tsv", "w") as out_fh:
	for i in range(0, len(probs_pleios)):
		print("{}\t{}".format(total_comb_pleios[i], probs_pleios[i]), file=out_fh)
	print("{}\t{}".format(Rest_comb, Prob_none), file=out_fh)
