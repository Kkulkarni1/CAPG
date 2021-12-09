#!/bin/bash

TEST_ONLY=0
EXE=roshan
NCPU=7
TARGET_FILE=/path/targets.txt
CMD=/path/ksd-roshan_wgs/roshan
GENOME_FILE=/path/combine.sam
FSA_FILES="/path/tet_A.fa /path/tet_B.fa"
SAMA_DIR=/path/subset_sam_A
SAMB_DIR=/path/subset_sam_B
OUTDIR=/path/err_files
VCFA=/path/indv_vcf_A
VCFB=/path/indv_vcf_B
INFOA=/path/info_A
INFOB=/path/info_B

function wait_on {
        NUM=`pgrep -c $EXE`
        while [ $NUM -ge $NCPU ]
        do
                echo "$NUM $EXE running..."
                sleep 2
                NUM=`pgrep -c $EXE`
        done
}

for GENOTYPE in SRR4124062 SRR4124066 SRR4124068 SRR4124074 SRR4124078 SRR8361734 SRR8361735 SRR8361736 SRR8361737 SRR8361738 SRR8736998 SRR8737008 SRR8737061 SRR8737062;
do
	for TARGET in `seq 1 1000`;
	do
		REF_NAMES=`grep "^\$TARGET\b" $TARGET_FILE | awk '{print \$2, \$3}'`
		RUN_THIS="$CMD --geno $GENOME_FILE --fsa_files $FSA_FILES --sam_files $SAMA_DIR/${GENOTYPE}_A.subset.sam $SAMB_DIR/${GENOTYPE}_B.subset.sam --ref_names $REF_NAMES --j extracted.fsa --vcf_files $VCFA/${GENOTYPE}_${TARGET}-A.vcf $VCFB/${GENOTYPE}_${TARGET}-B.vcf --o $INFOA/${GENOTYPE}-A.info.txt $INFOB/${GENOTYPE}-B.info.txt --name ${GENOTYPE} --soft-clipped 5 --po"
		if [ $TEST_ONLY -ge 1 ]
		then
			if [ ! -e $OUTDIR/${GENOTYPE} ]
			then
				echo "mkdir --parents $OUTDIR/${GENOTYPE}/"
			fi
			echo "$RUN_THIS 2> $OUTDIR/${GENOTYPE}/${GENOTYPE}_${TARGET}.err"
		else
			if [ ! -e $OUTDIR/${GENOTYPE} ]
			then
				mkdir --parents $OUTDIR/${GENOTYPE}/
			fi
			echo "$RUN_THIS 2> $OUTDIR/${GENOTYPE}/${GENOTYPE}_${TARGET}.err"
			wait_on
			$RUN_THIS 2> $OUTDIR/${GENOTYPE}/${GENOTYPE}_${TARGET}.err &
		fi
	done
done	 
