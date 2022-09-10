#!/bin/bash
# @author R. Kulkarni
# @file combine_infos.sh
#
# Purpose: Merge info files one per accession/region into a single info file per accession.

TEST_ONLY=0
CAPG_HOME=.
INDIR=${CAPG_HOME}/data/info
OUTDIR=${CAPG_HOME}/data/info

for GENOTYPE in SRR4124062	# this is a test script, using one accession
# SRR4124066 SRR4124068 SRR4124074 SRR4124078 SRR8361734 SRR8361735 SRR8361736 SRR8361737 SRR8361738 SRR8736998 SRR8737008 SRR8737061 SRR8737062;
do
	echo $GENOTYPE
	if [ $TEST_ONLY -ge 1 ]
	then
		echo "awk -F '\t' 'FNR>1 || NR==1' $INDIR/${GENOTYPE}_*_info.txt > $OUTDIR/${GENOTYPE}_info.txt"
	else
		awk -F '\t' 'FNR>1 || NR==1' $INDIR/${GENOTYPE}_*_info.txt > $OUTDIR/${GENOTYPE}_info.txt
	fi
done
