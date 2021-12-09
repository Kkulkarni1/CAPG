#!/usr/bin/python3
#
# Purpose: Extract truth state of each position to verify truth files truth_het_*.txt, truth_homeo_*.txt and truth_homoA_B_*.txt.
#
# Information: We call various truths, two at the individual level and two at the sample (multiple genotypes sampled from a population) level:
# * heterozygous loci in A and B subgenomes ==> truth_het_*.txt
# * homoeologous site, e.g. AA/CC, within individual; note this may not actually be homoeologous, as it could be segregating in either genome ==> truth_het_*.txt
# * true homoeologous site across individuals, e.g., all individuals are AA/CC ==> truth_homeo_*.txt
# * homologous site in A and B subgenomes across individuals (better called invariant), e.g. all individuals are AA/* in subgenome A ==> truth_homoA_B_*.txt
#
# Output files: truth_het_*.txt, truth_homeo_*.txt, truth_homoA_B_*.txt

do = "all"	#"homoeos"|"homos"|"heteros"		# limit was is done: all for do everything
write_files = False	# write truth files

from Bio import SeqIO
import pandas as pd
import numpy as np
import sys

# simulation conditions: hard-coded
homoeo_rates = ['0.005', '0.007', '0.010']
num_indiv = 50
genome_len = 100000

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

for hr in homoeo_rates:
	dir = "homr" + hr
	fsas = []
	for i in range(num_indiv):
		fsafile = dir + "/fasta_files/aln" + str(i) + ".fsa"
		seqs = []
		for sr in SeqIO.parse(fsafile, "fasta"):
			seqs.append(sr.seq)
		fsas.append(seqs)
	if do == "all" or do == "heteros":
		if write_files:
			fp_heteros = open("truth_het_" + hr + ".txt", "w")
			print("Genotype\tPosition\tTruth_het_A\tTruth_het_B\tTruth_homeo", file = fp_heteros)
		df_heteros = pd.DataFrame(columns = ['Genotype', 'Position', 'Truth_het_A', 'Truth_het_B', 'Truth_homeo'])
	if do == "all" or do == "homoeos":
		if write_files:
			fp_homoeos = open("truth_homeo_" + hr + ".txt", "w")
			print("Position\tTruth_homeo", file = fp_homoeos)
		df_homoeos = pd.DataFrame(columns = ['Position', 'Truth_homeo'])
	if do == "all" or do == "homos":
		if write_files:
			fp_homos = open("truth_homoA_B_" + hr + ".txt", "w")
			print("Position\tTruth_A\tTruth_B", file = fp_homos)
		df_homos = pd.DataFrame(columns = ['Position', 'Truth_A', 'Truth_B'])
	for j in range(genome_len):
		are_homoeo = 1
		are_homoA = 0
		are_homoB = 0
		for i in range(num_indiv):
			is_heteroA = 0
			is_heteroB = 0
			is_homoeo = 0
			seqs = fsas[i]
			if are_homoeo and (do == "all" or do == "homoeos" or (seqs[0][j] == seqs[2][j] or seqs[1][j] == seqs[3][j] or seqs[0][j] != seqs[1][j] or seqs[2][j] != seqs[3][j])):
				are_homoeo = 0
				if do == "homoeos":
					break
			if do == "all" or do == "homoeos" or do == "heteros":
				if seqs[0][j] != seqs[2][j] and seqs[0][j] == seqs[1][j] and seqs[2][j] == seqs[3][j]:
					is_homoeo = 1
				if seqs[0][j] != seqs[1][j]:
					is_heteroA = 1
				if seqs[2][j] != seqs[3][j]:
					is_heteroB = 1
			if  (do == "all" or do == "homos") and (are_homoA == 0 or are_homoB == 0):
				if i == 0:
					gen_A = seqs[0][j] + seqs[1][j]
					gen_B = seqs[2][j] + seqs[3][j]
				else:
					if are_homoA == 0 and seqs[0][j] + seqs[1][j] != gen_A:
						are_homoA = 1
						if do == "homos":
							break
					if are_homoB == 0 and seqs[2][j] + seqs[3][j] != gen_B:
						are_homoB = 1
						if do == "homos":
							break
			if do == "all" or do == "heteros":
				if write_files:
					print("aln" + str(i) + "\t" + str(j+1) + "\t" + str(is_heteroA) + "\t" + str(is_heteroB) + "\t" + str(is_homoeo), file = fp_heteros)
				df_heteros = df_heteros.append({'Genotype':'aln' + str(i), 'Position' : int(j+1), 'Truth_het_A' : is_heteroA, 'Truth_het_B' : is_heteroB, 'Truth_homeo' : is_homoeo}, ignore_index = True)
		if do == "all" or do == "homoeos":
			if write_files:
				print(str(j+1) + "\t" + str(are_homoeo), file = fp_homoeos)
			df_homoeos = df_homoeos.append({'Position' : int(j+1), 'Truth_homeo' : are_homoeo}, ignore_index = True)
		if do == "all" or do == "homos":
			if write_files:
				print(str(j+1) + "\t" + str(are_homoA) + "\t" + str(are_homoB), file = fp_homos)
			df_homos = df_homos.append({'Position' : int(j+1), 'Truth_A' : are_homoA, 'Truth_B' : are_homoB}, ignore_index = True)
	if do == "all" or do == "heteros":
		df_heteros.sort_values(by = ['Position','Genotype'], inplace = True, ignore_index = True)
		df_heteros_orig = pd.read_csv("truth_het_" + hr + ".txt", sep = "\t", dtype = {'Position' : np.int32})
		df_heteros_orig.sort_values(by = ['Position','Genotype'], inplace = True, ignore_index = True)
		print(df_heteros_orig.compare(df_heteros, keep_shape = False, keep_equal = True))
		print("Number of differences in A subgenome:" + str(sum(df_heteros['Truth_het_A'] - df_heteros_orig['Truth_het_A'])))
		print("Number of differences in B subgenome:" + str(sum(df_heteros['Truth_het_B'] - df_heteros_orig['Truth_het_B'])))
		print("Number of differences in AA/CC calls:" + str(sum(df_heteros['Truth_homeo'] - df_heteros_orig['Truth_homeo'])))
	if do == "all" or do == "homoeos":
		df_homoeos.sort_values(by = ['Position'], inplace = True, ignore_index = True)
		df_homoeos_orig = pd.read_csv("truth_homeo_" + hr + ".txt", sep = "\t")
		df_homoeos_orig.drop(0, inplace = True)
		df_homoeos_orig.sort_values(by = ['Position'], inplace = True, ignore_index = True)
		print(df_homoeos_orig.compare(df_homoeos, keep_shape = False, keep_equal = True))
		print("Number of differences:" + str(sum(df_homoeos['Truth_homeo'] - df_homoeos_orig['Truth_homeo'])))
	if do == "all" or do == "homos":
		df_homos.sort_values(by = ['Position'], inplace = True, ignore_index = True)
		df_homos_orig = pd.read_csv("truth_homoA_B_" + hr + ".txt", sep = "\t")
		df_homos_orig.drop(0, inplace = True)
		df_homos_orig.sort_values(by = ['Position'], inplace = True, ignore_index = True)
		print(df_homos_orig.compare(df_homos, keep_shape = False, keep_equal = True))
		print("Number of differences in A subgenome:" + str(sum(df_homos['Truth_A'] - df_homos_orig['Truth_A'])))
		print("Number of differences in B subgenome:" + str(sum(df_homos['Truth_B'] - df_homos_orig['Truth_B'])))
	if write_files:
		if do == "all" or do == "heteros":
			fp_heteros.close()
		if do == "all" or do == "homoeos":
			fp_homoeos.close()
		if do == "all" or do == "homos":
			fp_homos.close()
