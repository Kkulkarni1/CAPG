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
    #data_merge = data_merge.fillna(0)
    #data_merge["PP_homeo"] = (data_merge["PP_homeo"] * -1)
    
    data_PP_homeo = data_merge[["Position", "Truth_homeo", "PP_homeo"]]
    print(data_PP_homeo)

    if args.out_file:
        data_PP_homeo.to_csv(OUTFILE, header=["Position", "Truth_homeo", "PP_homeo"], index=None, sep='\t', mode='w')
