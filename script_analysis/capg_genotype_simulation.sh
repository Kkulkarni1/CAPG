#!/usr/bin/bash
# Purpose: Perform genotyping on simulated reads using CAPG.
# INPUT: SAM files generated by simulation.sh
# OUTPUT: err files generated by CAPG, individual VCF files

TEST=0					# test without doing anything
NCPU=7					# number of CPUs to use
OVERWRITE=1				# overwrite existing files
CAPG_HOME=.				# CAPG home directory
CAPG=$CAPG_HOME/script_main/wgs/release/capg_wgs	# CAPG executable
SIM_DIR=$CAPG_HOME/data/simulation	# simulation data
REF_SAM=ref.sam				# INPUT: alignment of subgenomic references
REFA_FSA=refA.fsa			# INPUT: subgenome reference A FASTA file
REFB_FSA=refB.fsa			# INPUT: subgenome reference B FASTA file
SAM_DIR=sam				# INPUT: where are sam file inputs
ERR_DIR=err				# OUTPUT: where to put stderr captured output
VCF_DIR=vcf				# OUTPUT: where to put vcf files
EXTRACTED_FSA_DIR=extracted		# OUTPUT: where to temporarily put extracted reference files
INFO_DIR=info				# [NOT USED] OUTPUT: where to put info files
HR_RATES="0.005 0.007"			# SIMULATOR: homoeologous SNP rate (manuscript: "0.005 0.007 0.100")
COV_RATES="5 10"			# SIMULATOR: coverage rates per chromosome (manuscript: "5 10 20" for hr = 0.007)
MM_RATES="0.000 0.001 0.010"		# SIMULATOR: mismatch rate (manuscript w/ mistake: "0.00 0.00 0.01" labeled as "0.00 0.01 0.10")
NREP=1					# SIMULATOR: number of genotypes to simulate (50)


EXE=$(basename $CAPG)

function wait_on {
        NUM=`pgrep -c $EXE`
        while [ $NUM -ge $NCPU ]
        do
                echo "$NUM $EXE running..."
                sleep 2
                NUM=`pgrep -c $EXE`
        done
}

for hr in $HR_RATES
do
	for cov in $COV_RATES
	do
		for mm in $MM_RATES
		do
			SIM_BASE="$SIM_DIR/homr${hr}mm${mm}"
			SIM_OUT="$SIM_BASE/cov${cov}"

			if [ $TEST == "0" -a ! -d "$SIM_OUT/$VCF_DIR" ]; then
				mkdir --parents $SIM_OUT/$VCF_DIR
			elif [ ! -d "$SIM_OUT/$VCF_DIR" ]; then
				echo "mkdir --parents $SIM_OUT/$VCF_DIR"
			fi
			if [ $TEST == "0" -a ! -d "$SIM_OUT/$ERR_DIR" ]; then
				mkdir --parents $SIM_OUT/$ERR_DIR
			elif [ ! -d "$SIM_OUT/$ERR_DIR" ]; then
				echo "mkdir --parents $SIM_OUT/$ERR_DIR"
			fi
			if [ $TEST == "0" -a ! -d "$SIM_OUT/$EXTRACTED_FSA_DIR" ]; then
				mkdir --parents $SIM_OUT/$EXTRACTED_FSA_DIR
			elif [ ! -d "$SIM_OUT/$EXTRACTED_FSA_DIR" ]; then
				echo "mkdir --parents $SIM_OUT/$EXTRACTED_FSA_DIR"
			fi

			# extract reference names from reference file
			REF_NAME_A=`head -n 1 $SIM_BASE/$REF_SAM | cut -f 2 | sed 's/SN://'`
			REF_NAME_B=`tail -n +2 $SIM_BASE/$REF_SAM | cut -f 1`

			for NUM in `seq 1 $NREP`
			do
				NUM=$((NUM-1))
				GENOTYPE=aln$NUM
				CMD="$CAPG --geno $SIM_BASE/$REF_SAM --fsa_files $SIM_BASE/$REFA_FSA $SIM_BASE/$REFB_FSA --sam_files $SIM_OUT/$SAM_DIR/${GENOTYPE}A.sam $SIM_OUT/$SAM_DIR/${GENOTYPE}B.sam --ref_names $REF_NAME_A $REF_NAME_B -j $SIM_OUT/$EXTRACTED_FSA_DIR/${GENOTYPE}_extracted.fsa --vcf_files $SIM_OUT/$VCF_DIR/${GENOTYPE}A.vcf $SIM_OUT/$VCF_DIR/${GENOTYPE}B.vcf --name ${GENOTYPE} -po" #-o $INFO_A/${GENOTYPE}A.txt $INFO_B/${GENOTYPE}B.txt"
				if [ ! -f $VCFA/${GENOTYPE}A.vcf -o $OVERWRITE -ge 1 ]; then
					wait_on
					echo "$CMD 2> $SIM_OUT/$ERR_DIR/${GENOTYPE}.err"
					if [ $TEST = "0" ]; then
						$CMD 2> $SIM_OUT/$ERR_DIR/${GENOTYPE}.err &
					fi
				fi
			done
		done
	done
done

# just cycle until done
NCPU=1
wait_on
