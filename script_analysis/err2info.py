#!/usr/bin/python3

import sys
import re
import argparse
from decimal import Decimal

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--in_file", required=True, nargs='+', help="input std.err file")
ap.add_argument("-o", "--out_file", required=False, help="info text file")
ap.add_argument("-e", "--equal_coverage", required=False, help="apply equal coverage test threshold at 0.05 threshold", action="store_true")
ap.add_argument("-c", "--coverage", required=False, help="coverage threshold", default = 8, type=int)
args = ap.parse_args()

for infile in args.in_file:
    INFILE = open(infile, 'rt')
    lines = INFILE.readlines()

    POSA=[]
    POSB=[]
    GENOTYPE_CALL=[]
    CALL_AGENOME=[]
    CALL_BGENOME = []
    PP_00 = []
    PP_01 = []
    PP_02 = []
    PP_10 = []
    PP_11 = []
    PP_12 = []
    PP_20 = []
    PP_21 = []
    PP_22 = []
    PP_MAX = []
    PPA_0 = []
    PPA_1 = []
    PPA_2 = []
    PPB_0 = []
    PPB_1 = []
    PPB_2 = []
    COV_A = []
    COV_B = []
    MJA = []
    MIA = []
    ECT = []
    ETA = []
    COUNT = 0
    for index, line in enumerate(lines):
        if "Sam files:" in line:	# extract accession number?
            genotype_match = re.match(r'.*\/([\w-]*.)_A.subset.sam', line)
            genotype = genotype_match.group(1)
            print(genotype)
        if "Reference names:" in line:	# extract chromosome name (without :start-end)
            chrom_match = re.match(r'.*Reference names:.*(arahy.*).*:.*(arahy.*).*:', line)
            chrom_A = chrom_match.group(1)
            print(chrom_A)
            chrom_B = chrom_match.group(2)
            print(chrom_B)
        if re.match(r'^Genotype ', line) is not None:
            POS_match = re.match(r'^Genotype \(\s*(\d*).\s*(\d*)', line)
            posA = POS_match.group(1)
            posB = POS_match.group(2)
            POSA.append(posA)
            POSB.append(posB)
            genotype_call_match = re.match(r'^Genotype.*\:.\s*(\D*).\(', line)
            genotype_call = genotype_call_match.group(1)
            GENOTYPE_CALL.append(genotype_call)
            call_Agenome_match = re.match(r'^Genotype.*\:.\s*(\w*).\D', line)
            call_Agenome = call_Agenome_match.group(1)
            CALL_AGENOME.append(call_Agenome)
            call_Bgenome_match = re.match(r'Genotype.*\/(\w*)', line)
            call_Bgenome = call_Bgenome_match.group(1)
            CALL_BGENOME.append(call_Bgenome)
            #print(lines[index-9])
            PP_00_match = re.match(r'.*\:\s*(\d.*)', lines[index - 9])
            pp_00 = Decimal(PP_00_match.group(1))
            PP_00.append(pp_00)
            PP_01_match = re.match(r'.*\:\s*(\d.*)', lines[index - 8])
            pp_01 = Decimal(PP_01_match.group(1))
            PP_01.append(pp_01)
            PP_02_match = re.match(r'.*\:\s*(\d.*)', lines[index - 7])
            pp_02 = Decimal(PP_02_match.group(1))
            PP_02.append(pp_02)
            PP_10_match = re.match(r'.*\:\s*(\d.*)', lines[index - 6])
            pp_10 = Decimal(PP_10_match.group(1))
            PP_10.append(pp_10)
            PP_11_match = re.match(r'.*\:\s*(\d.*)', lines[index - 5])
            pp_11 = Decimal(PP_11_match.group(1))
            PP_11.append(pp_11)
            PP_12_match = re.match(r'.*\:\s*(\d.*)', lines[index - 4])
            pp_12 = Decimal(PP_12_match.group(1))
            PP_12.append(pp_12)
            PP_20_match = re.match(r'.*\:\s*(\d.*)', lines[index - 3])
            pp_20 = Decimal(PP_20_match.group(1))
            PP_20.append(pp_20)
            PP_21_match = re.match(r'.*\:\s*(\d.*)', lines[index - 2])
            pp_21 = Decimal(PP_21_match.group(1))
            PP_21.append(pp_21)
            PP_22_match = re.match(r'.*\:\s*(\d.*)', lines[index - 1])
            pp_22 = Decimal(PP_22_match.group(1))
            PP_22.append(pp_22)

            PP_MAX.append(max(pp_00, pp_01, pp_02, pp_10, pp_11, pp_12, pp_20, pp_21, pp_22))

            MA_match = re.match(r'.*\s0\)\s(.)(.)', lines[index - 9])
            mja = MA_match.group(1)
            mia = MA_match.group(2)
            #print(mja, mia)
            MJA.append(mja)
            MIA.append(mia)
            ppA_0 = pp_00 + pp_01 + pp_02
            PPA_0.append(ppA_0)
            ppA_1 = pp_10 + pp_11 + pp_12
            PPA_1.append(ppA_1)
            ppA_2 = pp_20 + pp_21 + pp_22
            PPA_2.append(ppA_2)

            ppB_0 = pp_00 + pp_10 + pp_20
            PPB_0.append(ppB_0)
            ppB_1 = pp_01 + pp_11 + pp_21
            PPB_1.append(ppB_1)
            ppB_2 = pp_02 + pp_12 + pp_22
            PPB_2.append(ppB_2)

            COUNT += 1

        # WARNING: these lines MAY NOT be followed by a Genotype line, so you need to check and remove the last entry if the site is not genotyped
        if "Expected counts genome A" in line:
            cov = re.match(r'.*Expected counts genome A: (.*)$', line).group(1)
            covs = [float(i) for i in cov.split()]
            COV_A.append(sum(covs))
        if "Expected counts genome B" in line:
            cov = re.match(r'.*Expected counts genome B: (.*)$', line).group(1)
            covs = [float(i) for i in cov.split()]
            COV_B.append(sum(covs))
        # WARNING: here is the fix to remove those sites with coverage information but not genotype
        if "skip site" in line:
            COV_A.pop(-1)
            COV_B.pop(-1)
# INFO [capg.c::test_equal_homolog_coverage( 313)]: Equal coverage test (subgenome 1): eta = 0.586129, gamma1 = 1.000000, gamma2 = 0.894839 vs. eta = 0.586582, gamma1 = 0.000000, gamma2 = 0.500000; lrt = 13.543246 (-303.323377 -310.095000); pval = 2.331288e-04
        if "test_equal_homolog" in line:
            words = line.split(' ')
            ECT.append(float(words[-1]))
            ETA.append(float(words[-15].rstrip(',')))

    if args.out_file:
       OUTFILE = open(args.out_file, 'w')
       OUTFILE.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                "Genotype", "ChromA", "ChromB", "PositionA", "PositionB", "Genotype_call", "Call_Agenome", "Call_Bgenome", "Major allele", "Minor allele", "PP(0,0)",
                "PP(0,1)", "PP(0,2)", "PP(1,0)", "PP(1,1)", "PP(1,2)", "PP(2,0)", "PP(2,1)", "PP(2,2)", "PA(0)", "PA(1)",
                "PA(2)", "PB(0)", "PB(1)", "PB(2)", "CovA", "CovB", "ECTA","ETAA","ECTB","ETAB"))
    ect_index = 0
    filtered_covg = 0
    filtered_ect = 0
    for i in range(0, len(POSA)):
        #if (COV_A[i] >= 8 and COV_B[i] >=8 and PP_MAX[i] >= 0.98):
        #if (COV_A[i] >= 8 and COV_B[i] >=8 and PP_MAX[i] >= 0.95):
        pass_ectA = 1
        pass_ectB = 1
        alleles_a = list(CALL_AGENOME[i])
        alleles_b = list(CALL_BGENOME[i])
        if alleles_a[0] != alleles_a[1]:
            if args.equal_coverage and ECT[ect_index] < 0.05:
                pass_ectA = 0
        if alleles_b[0] != alleles_b[1]:
            if args.equal_coverage and ECT[ect_index] < 0.05:
                pass_ectB = 0
        if COV_A[i] >= args.coverage and COV_B[i] >= args.coverage and pass_ectA and pass_ectB:
           OUTFILE.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%f\t%f" % (
           genotype, chrom_A, chrom_B, POSA[i], POSB[i], GENOTYPE_CALL[i], CALL_AGENOME[i], CALL_BGENOME[i], MJA[i], MIA[i], PP_00[i], PP_01[i], PP_02[i], PP_10[i], PP_11[i], PP_12[i], PP_20[i], PP_21[i], PP_22[i], PPA_0[i], PPA_1[i], PPA_2[i], PPB_0[i], PPB_1[i], PPB_2[i], COV_A[i], COV_B[i]))
           if alleles_a[0] != alleles_a[1]:
               OUTFILE.write("\t%e\t%e" % (ECT[ect_index], ETA[ect_index]))
               ect_index += 1
           else:
               OUTFILE.write("\t%s\t%s" % ("NA","NA"))
           if alleles_b[0] != alleles_b[1]:
               OUTFILE.write("\t%e\t%e" % (ECT[ect_index], ETA[ect_index]))
               ect_index += 1
           else:
               OUTFILE.write("\t%s\t%s" % ("NA","NA"))
           OUTFILE.write("\n")
        else:
           if COV_A[i] < args.coverage or COV_B[i] < args.coverage:
               filtered_covg += 1
           if not pass_ectA or not pass_ectB:
               filtered_ect += 1
           if alleles_a[0] != alleles_a[1]:
               ect_index += 1
           if alleles_b[0] != alleles_b[1]:
               ect_index += 1
    print("Filtered", filtered_covg, "loci for coverage;", filtered_ect, "for ECT out of", len(POSA), "loci.")
