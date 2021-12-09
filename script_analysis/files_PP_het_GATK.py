import pandas as pd
import numpy as np
import argparse
ap = argparse.ArgumentParser()
ap.add_argument("-i1", "--in_file1", required=True, nargs='+', help="info file")
    #print(data_het_B_subset)
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

    data_merge = pd.merge(data, data_truth, how = 'left', left_on = ['Genotype', 'PositionA'], right_on = ['Genotype', 'Position'])

    data_merge = data_merge.fillna(0)
    print(data_merge)
    
    data_PR_het = data_merge[["Genotype", "Position","Truth_het_A", "Truth_het_B","PP_hetA", "PP_hetB"]]

    if args.out_file:
        data_PR_het.to_csv(OUTFILE, header=['Genotype', 'Position', 'Truth_het_A', 'Truth_het_B','PP_hetA', 'PP_hetB'], index=None, sep='\t', mode='w', na_rep='0')

