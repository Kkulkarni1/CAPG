#!/usr/bin/bash
# @file genotype_WGS.sh
# @author R. Kulkarni
#
# Purpose: run CAPG on sample peanut target
#
# Requires: compiled capg_wgs and various files in data directory, including several stub files, pointing to big data files elsewhere

TEST_ONLY=0
NCPU=7							# number of CPUs to use
NTARGET=10						# up to 1000 targets
CAPG_ROOT=.						# path to repository
TARGET_FILE=${CAPG_ROOT}/data/peanut/targets.txt	# list of target regions to genotype
CMD=${CAPG_ROOT}/script_main/wgs/release/capg_wgs	# CAPG executable
GENOME_FILE=${CAPG_ROOT}/data/peanut/peanut_ref.sam	# alignments of homoeologous targets
FSA_FILES="${CAPG_ROOT}/data/peanut/tet_A.fa ${CAPG_ROOT}/data/peanut/tet_B.fa"
							# subgenomic references NOT ON GITHUB
SAM_DIR=${CAPG_ROOT}/data/peanut/sam			# input directory for sam files
ERR_DIR=${CAPG_ROOT}/data/peanut/err			# output directory for stderr
VCF_DIR=${CAPG_ROOT}/data/peanut/vcf			# output directory for vcf files
EXTRACTED_DIR=${CAPG_ROOT}/data/peanut/extracted	# output directory for extracted target references

EXE=${CMD##*/}

function wait_on {
        NUM=`pgrep -c $EXE`
        while [ $NUM -ge $NCPU ]
        do
                echo "$NUM $EXE running..."
                sleep 2
                NUM=`pgrep -c $EXE`
        done
}

if [ ! -e $ERR_DIR/${GENOTYPE} ]
then
	echo "mkdir --parents $ERR_DIR"
	if [ $TEST_ONLY -ge 1 ]
	then
		mkdir --parents $ERR_DIR
	fi
fi
if [ ! -e $VCF_DIR ]
then
	echo "mkdir --parents $VCF_DIR"
	if [ $TEST_ONLY -ge 1 ]
	then
		mkdir --parents $VCF_DIR
	fi
fi
if [ ! -e $EXTRACTED_DIR ]
then
	echo "mkdir --parents $EXTRACTED_DIR"
	if [ $TEST_ONLY -ge 1 ]
	then
		mkdir --parents $EXTRACTED_DIR
	fi
fi

for GENOTYPE in SRR4124062	# this demo is for just one accession
# SRR4124066 SRR4124068 SRR4124074 SRR4124078 SRR8361734 SRR8361735 SRR8361736 SRR8361737 SRR8361738 SRR8736998 SRR8737008 SRR8737061 SRR8737062;
do
	for TARGET in `seq 1 $NTARGET`;
	do
		REF_NAMES=`grep "^\$TARGET\b" $TARGET_FILE | awk '{print \$2, \$3}'`

		RUN_THIS="$CMD --geno $GENOME_FILE --fsa_files $FSA_FILES --sam_files $SAM_DIR/${GENOTYPE}_A.subset.sam $SAM_DIR/${GENOTYPE}_B.subset.sam --ref_names $REF_NAMES --j $EXTRACTED_DIR/extracted_${GENOTYPE}_${TARGET}.fsa --name ${GENOTYPE} --equal -0 --vcf_files $VCF_DIR/${GENOTYPE}_${TARGET}-A.vcf $VCF_DIR/${GENOTYPE}_${TARGET}-B.vcf"
		if [ $TEST_ONLY -ge 1 ]
		then
			echo "$RUN_THIS 2> $ERR_DIR/${GENOTYPE}_${TARGET}.err"
		else
			echo "$RUN_THIS 2> $ERR_DIR/${GENOTYPE}_${TARGET}.err"
			wait_on
			$RUN_THIS 2> $ERR_DIR/${GENOTYPE}_${TARGET}.err &
		fi
	done
done
