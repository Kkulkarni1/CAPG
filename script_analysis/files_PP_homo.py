#!/usr/bin/python3
#
# Purpose: Merge truth and metric files for calling allelic SNPs (homoA, homoB).
#
# Input:
#	--in_file1	Metric files: /home/peanut/peanut2/CAPG/WGS/simulation/analysis_CAPG/PL_*_Ho.txt
#	--in_file2	Truth files: /home/peanut/peanut2/CAPG/WGS/simulation/truth_homo_*.txt
# Output:
#	--out_file	PR data files: /home/peanut/peanut2/CAPG/WGS/simulation/analysis_CAPG/files_plot/homozygous/PR_*_homo.txt

import pandas as pd
import numpy as np
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("-i1", "--in_file1", required=True, nargs='+', help="info file_1")
ap.add_argument("-i2", "--in_file2", required=True, nargs='+', help="truth file")
ap.add_argument("-o", "--out_file", required=True, help="PR_curve file")
args = ap.parse_args()

if args.out_file:
    OUTFILE = open(args.out_file, 'w')
for in_file1 in args.in_file1:
    data = pd.read_csv(in_file1, header = 0, sep = '\t')
    print(data)
for in_file2 in args.in_file2:
    data_truth = pd.read_csv(in_file2, header = 0, sep = '\t')
    print(data_truth)

    data_merge = pd.merge(data, data_truth, how='left', left_on = ['PositionA'], right_on = ['Position'])
    #data_merge.replace([np.inf, -np.inf], np.nan, inplace=True)
    data_merge = data_merge.fillna(0)
    data_merge["PP_homoA"] = (data_merge["PP_homoA"] * -1)
    data_merge["PP_homoB"] = (data_merge["PP_homoB"] * -1)
    
    data_PP_homo = data_merge[["Position", "Truth_A", "Truth_B", "PP_homoA", "PP_homoB"]]
    print(data_PP_homo)

    if args.out_file:
        data_PP_homo.to_csv(OUTFILE, header=["Position", "Truth_homoA", "Truth_homoB", "PP_homoA", "PP_homoB"], index=None, sep='\t', mode='w')
