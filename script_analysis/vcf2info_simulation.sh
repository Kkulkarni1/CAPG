#!/bin/bash
# @file vcf2info.sh
# @author R. Kulkarni
#
# Purpose:  Extract information for genotype calling from GATK vcf files.
#

TEST_ONLY=0
OVERWRITE=1
CAPG_HOME=.			# top directory of CAPG repository
SIM_DIR=data/simulation		# IN/OUTPUT: simulation data input/output
VCF_DIR=vcf			# INPUT: where vcf files are
INFO_DIR=info			# OUTPUT: info files
EXTENSION=_GATK_info.txt	# OUTPUT: extension to use for GATK
EXE_DIR=script_analysis		# directory with executables
EXE=vcf2info.py
HR_RATES="0.005 0.007"		# SIMULATOR: homoeologous SNP rate (manuscript: "0.005 0.007 0.100")
COV_RATES="5 10"		# SIMULATOR: coverage rates per chromosome (manuscript: "5 10 20" for hr = 0.007)
MM_RATES="0.000 0.001 0.010"	# SIMULATOR: mismatch rate (manuscript w/ mistake: "0.00 0.00 0.01" labeled as "0.00 0.01 0.10")
NREP=1				# SIMULATOR: number of genotypes to simulate (50)


for hr in $HR_RATES
do
	for mm in $MM_RATES
	do
		for cvg in $COV_RATES
		do
			DIR="$CAPG_HOME/$SIM_DIR/homr${hr}mm${mm}/cov${cvg}"

			if [ ! -d "$DIR" ]; then
				echo "$DIR missing..."
				continue
			fi

			if [ ! -d "$DIR/GATK/$INFO_DIR" ]; then
				echo "mkdir --parents $DIR/GATK/$INFO_DIR"
				if [ $TEST_ONLY == 0 ]; then
					mkdir --parents $DIR/GATK/$INFO_DIR
				fi
			fi
			
			for i in `seq 1 $NREP`; do
				i=$((i-1))
				if [ -f "$DIR/GATK/$VCF_DIR/sim${i}A.vcf" -a -f "$DIR/GATK/$VCF_DIR/sim${i}B.vcf" ]; then
					CMD="python $EXE_DIR/$EXE -i1 $DIR/GATK/$VCF_DIR/sim${i}A.vcf -i2 $DIR/GATK/$VCF_DIR/sim${i}B.vcf -o $DIR/GATK/$INFO_DIR/aln${i}${EXTENSION}"
					echo $CMD
					if [ $TEST_ONLY == 0 ]; then
						$CMD
					fi
				else
					echo "$DIR/GATK/$VCF_DIR/sim${i}A.vcf or $DIR/GATK/$VCF_DIR/sim${i}B.vcf missing..."
				fi
			done
		done
	done
done
