#!/usr/bin/python3

import pandas as pd
import numpy as np
import requests, sys
from itertools import combinations
from itertools import product
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

Pleiotropies=["Agon_early_early", "Agon_early_late", "Agon_early_early",
"Antagon_early_early", "Antagon_early_late", "Agon_early_early"]
PATH1=sys.argv[1]
PATH2=sys.argv[2]
probs_pleios = {}
total_comb_pleios = {}
diseases_pairs = {}
totals_disease = {}

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
for subdir in os.listdir(PATH1):
	pair_APs = []
	pair_dis = subdir.split("_")
	#Now we go with the count of OBSERVED combined
	#First with the paired
	INPUT_COMB_PATH=sys.argv[1] + \
	pair_dis[0] + "_" + pair_dis[1] + "/pleios.txt"
	Comb_pair_val = count_number_of_pleios(INPUT_COMB_PATH)
	INPUT_Single_COMB_PATH = sys.argv[1] + pair_dis[0] + "_" + pair_dis[1] + "/Age_threeshold_10/SinglePleiosLoc.tsv"
	single_comb = count_number_of_pleios(INPUT_Single_COMB_PATH) - 1
	Total_comb = Comb_pair_val + single_comb
	#total_comb_pleios[subdir] = Total_comb
	total_comb_pleios[subdir] = Total_comb

	#for each one alone
for subdir in os.listdir(PATH2):
	#Paired SNP pleiotropies
	INPUT_PATH=sys.argv[2] + subdir
	total_pair_disease = count_number_of_pleios(INPUT_PATH + "/pleios.txt")
	#Single SNP pleiotropies
	total_single_disease = count_number_of_pleios(INPUT_PATH +
	"/Age_threeshold_10/SinglePleiosLoc.tsv") - 1
	total_disease = total_pair_disease + total_single_disease
	#totals_disease[subdir] = total_disease
	totals_disease[subdir] = total_disease

#TOTAL PROBS COMBINATIONS
#combinations_list = [{j: totals_disease[j] for j in i} for iº in combinations(totals_disease.items(), 2)]
combinations_list = list(map(dict, combinations(totals_disease.items(), 2)))
for element_dict in combinations_list:
	key_values = []
	for key in element_dict:
		key_values.append(key)
	probs_pleios[key_values[0] + "_" + key_values[1]] = \
	(element_dict[key_values[0]]/sum(list(totals_disease.values())) * \
	element_dict[key_values[1]]/sum(list(totals_disease.values())))

#Prob_none = 1 - sum(probs_pleios)
#Rest_comb = total_pleios - sum(total_comb_pleios)

#Finally, we print the ouput in a new folders
with open(PATH1 + "freqs_total.tsv", "w") as out_fh:
	print(sum(list(totals_disease.values())))
	print(sum(list(total_comb_pleios.values())))
	print("{}\t{}\t{}".format("Combined_diseases", "Total_observed", "Exp_prob"), file=out_fh)
	for key in total_comb_pleios:
		list2 = key.split("_")
		if list2[0] == list2[1]:
			print(key)
			print((totals_disease[list2[0]]/sum(list(totals_disease.values()))))
			print("{}\t{}\t{}".format(key, total_comb_pleios[key],
			(totals_disease[list2[0]]/sum(list(totals_disease.values()))) * (totals_disease[list2[0]]/sum(list(totals_disease.values())))), file=out_fh)
		elif key in probs_pleios:
			print("{}\t{}\t{}".format(key, total_comb_pleios[key], probs_pleios[key]), file=out_fh)
		else:
			list2 = key.split("_")
			rev_list2 = list2[::-1]
			key2 = "_".join(rev_list2)
			print("{}\t{}\t{}".format(key, total_comb_pleios[key], probs_pleios[key2]), file=out_fh)
	#print("{}\t{}\t{}".format("Rest", Rest_comb, Prob_none), file=out_fh)
