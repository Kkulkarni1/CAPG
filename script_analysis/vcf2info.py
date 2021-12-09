#!/usr/bin/python3
#
# Purpose: Convert GATK vcf files into information files.
# Input: vcf files
# Output: Info files (processed data files)
# Usage: bash ~/scripts/vcf2info.sh

import re
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("-i1", "--in_file_A", required=True, nargs='+', help="input vcf file A")
ap.add_argument("-i2", "--in_file_B", required=True, nargs='+', help="input vcf file B")
ap.add_argument("-o", "--out_file", required=False, help="info text file")
args = ap.parse_args()

if args.out_file:
    OUTFILE = open(args.out_file, 'w')
    #OUTFILE.write("%s\t%s\t%s\t%s\t%s\n" % ("Genotype", "Position", "PLA(0)", "PLA(1)", "PLA(2)"))
    OUTFILE.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("Genotype", "ChromA", "ChromB", "PositionA", "PositionB", "Genotype_Call", "Ref_A_allele", "Alt_A_allele", "Ref_B_allele", "Alt_B_allele", "PLA(0)", "PLA(1)", "PLA(2)", "PLB(0)", "PLB(1)", "PLB(2)", "CovA", "CovB"))

#INFILE = open("aln0-A.vcf",'rt')
for infile_A in args.in_file_A:
    INFILE_A = open(infile_A, 'rt')
    lines_A = INFILE_A.readlines()

for infile_B in args.in_file_B:
    INFILE_B = open(infile_B, 'rt')
    lines_B = INFILE_B.readlines()

for line_A, line_B in zip(lines_A, lines_B):
    line_A = line_A.strip('\n')
    line_B = line_B.strip('\n')
    GENOTYPE_CALL_A = []
    GENOTYPE_CALL_B = []
    REF_ALLELE_A = []
    ALT_ALLELE_A = []
    REF_ALLELE_B = []
    ALT_ALLELE_B = []
    CHROM_A = []
    POSITION_A = []
    PL_REF_A = []
    PL_HET_A = []
    PL_ALT_A = []
    COVG_A = []
    CHROM_B = []
    POSITION_B = []
    PL_REF_B = []
    PL_HET_B = []
    PL_ALT_B = []
    COVG_B = []
    if "#CHROM" in line_A:
        GENOTYPE_MATCH_A = re.match(r'.*(aln.*).*', line_A)
        GENOTYPE_A = GENOTYPE_MATCH_A.group(1)
        print(GENOTYPE_A)

    # Genome_B	9924	.	C	<NON_REF>	.	.	.	GT:AD:DP:GQ:PL	0/0:49,0:49:99:0,120,1800
    if "Genome_A" in line_A and "##" not in line_A and "SB" not in line_A and ":GQ" in line_A:
        CHROM_A.append(line_A.split('\t')[0])
        POSITION_A.append(line_A.split('\t')[1])
        GENOTYPE_CALL = re.split('[\t : ]', line_A)[13]
        GENOTYPE_CALL_A.append(sum([int(i) for i in re.split('[/|]', GENOTYPE_CALL)]))
        COVG_A.append(re.split('[\t : ]', line_A)[15])
        REF_ALLELE_A.append(re.split('[\t , ]', line_A)[3])
        ALT_ALLELE_A.append(re.split('[\t , ]', line_A)[4])
        raw_PL_A = re.split('[\t : ]', line_A)[17]
        #print(raw_PL)
        PL_REF_A.append(raw_PL_A.split(',')[0])
        PL_HET_A.append(raw_PL_A.split(',')[1])
        PL_ALT_A.append(raw_PL_A.split(',')[2])

    #Genome_B	832	.	T	C,<NON_REF>	260.64	.	BaseQRankSum=-0.345;DP=36;ExcessHet=3.0103;MLEAC=1,0;MLEAF=0.500,0.00;MQRankSum=0.000;RAW_MQandDP=129600,36;ReadPosRankSum=0.812	GT:AD:DP:GQ:PL:SB	0/1:26,9,0:35:99:268,0,917,346,944,1290:12,14,6,3
    elif "Genome_A" in line_A and "##" not in line_A and ":SB" in line_A and "PID" not in line_A and ":GQ" in line_A:
        #print(line)
        CHROM_A.append(line_A.split('\t')[0])
        POSITION_A.append(line_A.split()[1])
        GENOTYPE_CALL = re.split('[\t : ]', line_A)[14]
        GENOTYPE_CALL_A.append(sum([int(i) for i in re.split('[/|]', GENOTYPE_CALL)]))
        COVG_A.append(re.split('[\t : ]', line_A)[16])
        REF_ALLELE_A.append(re.split('[\t , ]', line_A)[3])
        ALT_ALLELE_A.append(re.split('[\t , ]', line_A)[4])
        raw_PL_A = line_A.split()[9]
        #raw_PL_1 = re.split('[,:]', raw_PL)[]
        PL_REF_A.append(re.split('[,:]', raw_PL_A)[6])
        PL_HET_A.append(re.split('[,:]', raw_PL_A)[7])
        PL_ALT_A.append(re.split('[,:]', raw_PL_A)[8])

    # Genome_B	99801	.	G	C,<NON_REF>	134.64	.	BaseQRankSum=-1.319;DP=21;ExcessHet=3.0103;MLEAC=1,0;MLEAF=0.500,0.00;MQRankSum=0.000;RAW_MQandDP=75600,21;ReadPosRankSum=-3.086	GT:AD:DP:GQ:PGT:PID:PL:PS:SB	0|1:15,6,0:21:99:0|1:99801_G_C:142,0,555,187,573,759:99801:1,14,0,6
    elif "Genome_A" in line_A and "##" not in line_A and ":SB" in line_A and "PID" in line_A and ":GQ" in line_A:
        CHROM_A.append(line_A.split('\t')[0])
        POSITION_A.append(line_A.split()[1])
        GENOTYPE_CALL = re.split('[\t : ]', line_A)[17]
        GENOTYPE_CALL_A.append(sum([int(i) for i in re.split('[/|]', GENOTYPE_CALL)]))
        COVG_A.append(re.split('[\t : ]', line_A)[19])
        REF_ALLELE_A.append(re.split('[\t , ]', line_A)[3])
        ALT_ALLELE_A.append(re.split('[\t , ]', line_A)[4])
        raw_PL_A = line_A.split()[9]
        PL_REF_A.append(re.split('[,:]', raw_PL_A)[8])
        PL_HET_A.append(re.split('[,:]', raw_PL_A)[9])
        PL_ALT_A.append(re.split('[,:]', raw_PL_A)[10])
    # Genome_B 42046 . C G,<NON_REF> 0.01 . DP=1;MLEAC=0,0;MLEAF=NaN,NaN;RAW_MQandDP=3600,1 GT ./.
    elif "Genome_A" in line_A and "##" not in line_A and "DP=" in line_A:
        CHROM_A.append(line_A.split('\t')[0])
        POSITION_A.append(line_A.split()[1])
        GENOTYPE_CALL_A.append('0')
        COVG_A.append(re.split('[\t ; = ]', line_A)[8])
        REF_ALLELE_A.append(re.split('[\t , ]', line_A)[3])
        ALT_ALLELE_A.append(re.split('[\t , ]', line_A)[4])
        PL_REF_A.append(0)
        PL_HET_A.append(0)
        PL_ALT_A.append(0)
    # Genome_B	12777	.	T	A,<NON_REF>	0.01	.	MLEAC=0,0;MLEAF=NaN,NaN	GT:PGT:PID:PS	.|.:0|1:12777_T_A:12777
    elif "Genome_A" in line_A and "##" not in line_A and "MLEAF=NaN" in line_A:
        CHROM_A.append(line_A.split('\t')[0])
        POSITION_A.append(line_A.split()[1])
        GENOTYPE_CALL_A.append('0')
        COVG_A.append(0)
        REF_ALLELE_A.append(re.split('[\t , ]', line_A)[3])
        ALT_ALLELE_A.append(re.split('[\t , ]', line_A)[4])
        PL_REF_A.append(0)
        PL_HET_A.append(0)
        PL_ALT_A.append(0)

    if "#CHROM" in line_B:
        GENOTYPE_MATCH_B = re.match(r'.*(aln.*).*', line_B)
        GENOTYPE_B = GENOTYPE_MATCH_B.group(1)
        print(GENOTYPE_B)
    if "Genome_B" in line_B and "##" not in line_B and "SB" not in line_B and ":GQ" in line_B:
        #print(line)
        CHROM_B.append(line_B.split('\t')[0])
        POSITION_B.append(line_B.split('\t')[1])
        REF_ALLELE_B.append(re.split('[\t , ]', line_B)[3])
        ALT_ALLELE_B.append(re.split('[\t , ]', line_B)[4])
        GENOTYPE_CALL = re.split('[\t : ]', line_B)[13]
        GENOTYPE_CALL_B.append(sum([int(i) for i in re.split('[/|]', GENOTYPE_CALL)]))
        COVG_B.append(re.split('[\t : ]', line_B)[15])
        raw_PL_B = re.split('[\t : ]', line_B)[17]
        #print(raw_PL)
        PL_REF_B.append(raw_PL_B.split(',')[0])
        PL_HET_B.append(raw_PL_B.split(',')[1])
        PL_ALT_B.append(raw_PL_B.split(',')[2])
    elif "Genome_B" in line_B and "##" not in line_B and ":SB" in line_B and "PID" not in line_B and ":GQ" in line_B:
        #print(line)
        CHROM_B.append(line_B.split('\t')[0])
        POSITION_B.append(line_B.split()[1])
        REF_ALLELE_B.append(re.split('[\t , ]', line_B)[3])
        ALT_ALLELE_B.append(re.split('[\t , ]', line_B)[4])
        GENOTYPE_CALL = re.split('[\t : ]', line_B)[14]
        GENOTYPE_CALL_B.append(sum([int(i) for i in re.split('[/|]', GENOTYPE_CALL)]))
        COVG_B.append(re.split('[\t : ]', line_B)[16])
        raw_PL_B = line_B.split()[9]
        #raw_PL_1 = re.split('[,:]', raw_PL)[]
        PL_REF_B.append(re.split('[,:]', raw_PL_B)[6])
        PL_HET_B.append(re.split('[,:]', raw_PL_B)[7])
        PL_ALT_B.append(re.split('[,:]', raw_PL_B)[8])
    elif "Genome_B" in line_B and "##" not in line_B and ":SB" in line_B and "PID" in line_B and ":GQ" in line_B:
        CHROM_B.append(line_B.split('\t')[0])
        POSITION_B.append(line_B.split()[1])
        REF_ALLELE_B.append(re.split('[\t , ]', line_B)[3])
        ALT_ALLELE_B.append(re.split('[\t , ]', line_B)[4])
        GENOTYPE_CALL = re.split('[\t : ]', line_B)[17]
        GENOTYPE_CALL_B.append(sum([int(i) for i in re.split('[/|]', GENOTYPE_CALL)]))
        COVG_B.append(re.split('[\t : ]', line_B)[19])
        raw_PL_B = line_B.split()[9]
        PL_REF_B.append(re.split('[,:]', raw_PL_B)[8])
        PL_HET_B.append(re.split('[,:]', raw_PL_B)[9])
        PL_ALT_B.append(re.split('[,:]', raw_PL_B)[10])
    elif "Genome_B" in line_B and "##" not in line_B and "DP=" in line_B:
        CHROM_B.append(line_B.split('\t')[0])
        POSITION_B.append(line_B.split()[1])
        GENOTYPE_CALL_A.append('0')
        REF_ALLELE_B.append(re.split('[\t , ]', line_B)[3])
        ALT_ALLELE_B.append(re.split('[\t , ]', line_B)[4])
        COVG_B.append(re.split('[\t ; = ]', line_B)[8])
        PL_REF_B.append(0)
        PL_HET_B.append(0)
        PL_ALT_B.append(0)
    elif "Genome_B" in line_B and "##" not in line_B and "MLEAF=NaN" in line_B:
        CHROM_B.append(line_B.split('\t')[0])
        POSITION_B.append(line_B.split()[1])
        GENOTYPE_CALL_A.append('0')
        REF_ALLELE_B.append(re.split('[\t , ]', line_B)[3])
        ALT_ALLELE_B.append(re.split('[\t , ]', line_B)[4])
        COVG_B.append(re.split('[\t ; = ]', line_B)[8])
        PL_REF_B.append(0)
        PL_HET_B.append(0)
        PL_ALT_B.append(0)

    if len(POSITION_A) != len(POSITION_B):
        sys.exit("Mismatched file lengths.\n")


    if args.out_file:
        for i in range(0, len(POSITION_A)):
            #print(i, GENOTYPE_A, POSITION_A[i], POSITION_B[i], GENOTYPE_CALL[i], REF_ALLELE[i], ALT_ALLELE[i], PL_REF_A[i], PL_HET_A[i], PL_ALT_A[i], PL_REF_B[i], PL_HET_B[i], PL_ALT_B[i], COVG_A[i], COVG_B[i])
            if ALT_ALLELE_A[i] == "<NON_REF>":
                ALT_ALLELE_A[i] = 'NA'
            if ALT_ALLELE_B[i] == "<NON_REF>":
                ALT_ALLELE_B[i] = 'NA'
            OUTFILE.write("%s\t%s\t%s\t%s\t%s\t%s/%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (GENOTYPE_A, CHROM_A[i], CHROM_B[i], POSITION_A[i], POSITION_B[i], GENOTYPE_CALL_A[i], GENOTYPE_CALL_B[i], REF_ALLELE_A[i], ALT_ALLELE_A[i], REF_ALLELE_B[i], ALT_ALLELE_B[i], PL_REF_A[i], PL_HET_A[i], PL_ALT_A[i], PL_REF_B[i], PL_HET_B[i], PL_ALT_B[i], COVG_A[i], COVG_B[i]))
