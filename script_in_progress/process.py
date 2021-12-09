import os
import sys
import argparse
import pandas as pd
import dataprocessing as dp

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--in_file", required=True, nargs="+", help="input info file")
ap.add_argument("-o", "--out_file", required=True, nargs=1, help="output name root")
ap.add_argument(
    "-a",
    "--a_ref",
    required=False,
    nargs=1,
    help='print "A" genome vcf using supplied reference fasta',
)
ap.add_argument(
    "-b",
    "--b_ref",
    required=False,
    nargs=1,
    help='print "B" genome vcf using supplied reference fasta',
)
ap.add_argument(
    "-p",
    "--pseudo",
    action="store_true",
    help="print pseudo-vcf based on A+B reference fasta",
)

# ap.add_argument('-l', '--log', help='log_10 transform PP values', action='store_true')
args = ap.parse_args()
outfile = args.out_file[0]

if args.pseudo and not (args.a_ref and args.b_ref):
    sys.exit("error: %s\n" % ("Must use --a_ref and --b_ref with pseudo"))
vcf_args = [False, False, False, False]
if args.a_ref or args.b_ref:
    vcf_args = [args.out_file, args.a_ref, args.b_ref, args.pseudo]

f = dp.Readfiles(args.in_file)
#setup het* file printer
phet = dp.Het(outfile)

#load files
ds = f.topandas()
#calculate pp
cpp = dp.pp(ds)
ds = cpp.calc()
ds = (
    ds[
        [
            "ChromA",
            "ChromB",
            "PositionA",
            "PositionB",
            "Genotype",
            "Genotype_call",
            "Call_Agenome",
            "Call_Bgenome",
            "pp_homeo",
            "pp_homo_a",
            "pp_homo_b",
            "pp_het_a",
            "pp_het_b",
            "CovA",
            "CovB",
        ]
    ]
    .sort_values(by=["ChromA", "ChromB", "PositionA", "PositionB", "Genotype"])
    .rename(
        {
            "pp_homeo": "PP(homeo)",
            "pp_homo_a": "PP(homoA)",
            "pp_homo_b": "PP(homoB)",
            "pp_het_a": "PP(hetA)",
            "pp_het_b": "PP(hetB)",
        },
        axis=1,
    )
)

#print het files
phet.printHetFiles(ds)
#print vcf
if args.a_ref or args.b_ref:
    active_chrom = ds["ChromA"].unique().tolist()
    active_chrom += ds["ChromB"].unique().tolist()
    vcf_args.append(active_chrom)
    v = dp.Vcf(vcf_args)
    v.print_vcf_from_dataframe(ds)

