#!/bin/bash
# @author Yudi Zhang, Karin S. Dorman, Roshan Kulkarni
# @file simulate_data.sh
#
# Purpose: Simulate WGS data. This does 5% of the simulation run in CAPG
# manuscript but in a full factorial design over homoeologous rates, coverage
# rates, and subgenome reference mismatche rate.
#

TEST=0					# test this code without doing anything
NCPU=7					# number of CPU allocated to this script

OVERWRITE_SUBGENOMES=1			# simulate subgenomes again
OVERWRITE_GENOMES=$OVERWRITE_SUBGENOMES	# simulate individual genomes again [WARNING: if overwriting subgenomes, it WILL overwrite sample!!]
OVERWRITE_READS=1			# simulate reads from subgenomes again
OVERWRITE_REFERENCES=1			# simulate references from subgenomes again
OVERWRITE_ALIGNMENTS=1			# simulate alignments from references and reads again

CAPG_HOME=.				# where input and output should go
SIMULATOR=$CAPG_HOME/src/wgs/capg_sim	# simulator executable		
SIM_DIR=$CAPG_HOME/data/simulation	# directory for simulation data
REFERENCE="$SIM_DIR/refA.fsa"		# SIMULATOR: FASTA with simulated reference A
ERR_FILE1="$SIM_DIR/miseq250R1.txt"	# ART: art_illumina error files	
ERR_FILE2="$SIM_DIR/miseq250R2.txt"
					# SIMULATOR: implement full factorial design on these variables
HR_RATES="0.005" # 0.007"	 		# SIMULATOR: homoeologous SNP rate (manuscript: "0.005 0.007 0.100")
COV_RATES="5" # 10"			# SIMULATOR: coverage rates per chromosome, twice this per subgenome (manuscript: "5 10 20" for hr = 0.007)
MM_RATES="0.000 0.001 0.010"		# SIMULATOR: mismatch rate (manuscript w/ mistake: "0.00 0.00 0.01" labeled as "0.00 0.01 0.10")
NREP=50					# SIMULATOR: number of genotypes to simulate (50)

# simulation conditions: fixed parameters
ALPHA=100		# SIMULATOR: allele proportions from Beta(ALPHA, BETA)
BETA=100
ALLELIC_RATE=0.005	# SIMULATOR: probability of allelic SNP at sites that are not homoeologous loci
READ_LEN=150		# ART: read length
FRAG_LEN_MEAN=300	# ART: fragment length mean
FRAG_LEN_SD=10		# ART: fragment length standard deviation

EXE=$(basename $SIMULATOR)

function wait_on_cpu {
	NUM=`pgrep -c $EXE`
	while [ $NUM -ge $NCPU ]
	do
		echo "$NUM $EXE running..."
		sleep 2
		NUM=`pgrep -c $EXE`
	done
}

function wait_on_finish {
	NUM=`pgrep -c $EXE`
	while [ $NUM -ge 1 ]
	do
		echo "$NUM $EXE running..."
		sleep 2
		NUM=`pgrep -c $EXE`
	done
}

for hr in $HR_RATES
do
	PARENT_DIR="$SIM_DIR/homr${hr}mm0.000"	# folder where true subgenomes and individual genomes simulated
	if [ $TEST == "0" -a ! -d "$PARENT_DIR/sample" ]; then
		mkdir --parents $PARENT_DIR/sample
	elif [ ! -d "$PARENT_DIR/sample" ]; then
		echo "mkdir --parents $PARENT_DIR/sample"
	fi

	# simulate subgenomes with given homoeologous SNP rate
	# and simulate individual genomes with given allelic SNP rate
	if [ ! -f "$PARENT_DIR/ref.fsa" -o $OVERWRITE_SUBGENOMES -ge 1 ]; then
		CMD="$SIMULATOR -e $REFERENCE -j $hr -g $ALLELIC_RATE -p 0.00 -a $ALPHA -b $BETA -s $RANDOM -o $PARENT_DIR/sample/indiv -f $PARENT_DIR/ref -r $PARENT_DIR/ref.sam -n $NREP"
        	# Simulating reference B genome, allotetraploid individuals and generating truth files
		echo "$CMD 2> $PARENT_DIR/truth_${hr}.txt"
		if [ $TEST == "0" ]; then
			$CMD 2> $PARENT_DIR/truth_${hr}.txt
		fi
	# just simulate individual genomes with given allelic SNP rate, taking subgenomes from existing reference file
	elif [ $OVERWRITE_GENOMES -ge 1 ]; then
		CMD="$SIMULATOR -e $PARENT_DIR/ref.fsa -j $hr -g $ALLELIC_RATE -p 0.00 -a $ALPHA -b $BETA -s $RANDOM -o $PARENT_DIR/sample/indiv -f $PARENT_DIR/ref -r $PARENT_DIR/ref.sam -n $NREP"
        	# Simulating reference B genome, allotetraploid individuals and generating truth files
		echo "$CMD 2> $PARENT_DIR/truth_${hr}.txt"
		if [ $TEST = "0" ]; then
			$CMD 2> $PARENT_DIR/truth_${hr}.txt
		fi
	fi

	# simulate reads
	EXE=art_illumina
	for cov in $COV_RATES
	do
		READ_DIR="$PARENT_DIR/cov${cov}/fastq"
		if [ $TEST == "0" -a ! -d "$READ_DIR" ]; then 
			mkdir --parents $READ_DIR
		elif [ ! -d "$READ_DIR/fastq" ]; then
			echo "mkdir --parents $READ_DIR"
		fi
		for m in `seq 1 $NREP` #number of individuals
		do
			m=$((m-1))
			# simulate reads
			if [ ! -f "$PARENT_DIR/fastq/sim${m}.1.fq -o ! -f "$PARENT_DIR/fastq/sim${m}.2.fq -o $OVERWRITE_READS -ge 1 ]; then
				echo "art_illumina -1 $ERR_FILE1 -2 $ERR_FILE2 --len $READ_LEN --in $PARENT_DIR/sample/indiv${m}.fsa --out $READ_DIR/sim${m}. --fcov ${cov} --paired --mflen $FRAG_LEN_MEAN --sdev $FRAG_LEN_SD --rndSeed $RANDOM"
				if [ $TEST == "0" ]; then
					wait_on_cpu
					art_illumina -1 $ERR_FILE1 -2 $ERR_FILE2 --len $READ_LEN --in $PARENT_DIR/sample/indiv${m}.fsa --out $READ_DIR/sim${m}. --fcov ${cov} --paired --mflen $FRAG_LEN_MEAN --sdev $FRAG_LEN_SD --rndSeed $RANDOM &
				fi
			fi
		done
		wait_on_finish
	done

	# simulate references and complete alignment of each set of reads to each reference
	for mm in $MM_RATES
	do
		REF_DIR="$SIM_DIR/homr${hr}mm${mm}"	# folder with references

		if [ $(echo "$mm > 0" | bc) -gt 0 ]; then
			if [ $TEST == "0" -a ! -d "$REF_DIR" ]; then
				mkdir --parents $REF_DIR
			elif [ ! -d "$REF_DIR" ]; then
				echo "mkdir --parents $REF_DIR"
			fi

			# simulate imperfect references with homoeologous and allelic SNPs and given homoeologous rate
			if [ ! -f "$REF_DIR/ref.fsa" -o $OVERWRITE_REFERENCES -ge 1 ]; then
				CMD="$SIMULATOR -e $PARENT_DIR/ref.fsa -j $hr -g $ALLELIC_RATE -p $mm -s $RANDOM -f $REF_DIR/ref -r $REF_DIR/ref.sam -n 0"
				# simulate reference A, B genome
				echo "$CMD 2> $REF_DIR/truth_${hr}.txt"
				if [ $TEST = "0" ]; then
					$CMD 2> $REF_DIR/truth_${hr}.txt
				fi
			fi
		fi

		# generate indices of the references
		if [ $TEST == "0" -a ! -f "$REF_DIR/refA.fsa.bwt" ]; then
			bwa index $REF_DIR/refA.fsa
		elif [ ! -f "$REF_DIR/refA.fsa.bwt" ]; then
			echo "bwa index $REF_DIR/refA.fsa"
		fi
		if [ $TEST == "0" -a ! -f "$REF_DIR/refB.fsa.bwt" ]; then
			bwa index $REF_DIR/refB.fsa
		elif [ ! -f "$REF_DIR/refB.fsa.bwt" ]; then
			echo "bwa index $REF_DIR/refB.fsa"
		fi

		# generate alignments to each (flawed) reference for each set of reads
		#EXE="bwa"
		for cov in $COV_RATES	# read coverage level
		do
			READ_DIR="$PARENT_DIR/cov${cov}/fastq"
			SAM_DIR="$REF_DIR/cov${cov}/sam"
			if [ $TEST == "0" -a ! -d "$SAM_DIR" ]; then
				mkdir --parents $SAM_DIR
			elif [ ! -d "$SAM_DIR" ]; then
				echo "mkdir --parents $SAM_DIR"
			fi

			for m in `seq 1 $NREP` # individuals
			do
				m=$((m - 1))
				# align reads to A and B subgenome references

				if [ ! -f "$SAM_DIR/aln${m}A.sam" -o $OVERWRITE_ALIGNMENTS -ge 1 ]; then
					echo "bwa mem -t $NCPU $REF_DIR/refA.fsa $READ_DIR/sim${m}.1.fq $READ_DIR/sim${m}.2.fq > $SAM_DIR/aln${m}A.sam"
					if [ $TEST == "0"  ]; then
						#wait_on_cpu
						bwa mem -t $NCPU $REF_DIR/refA.fsa $READ_DIR/sim${m}.1.fq $READ_DIR/sim${m}.2.fq > $SAM_DIR/aln${m}A.sam
					fi
				fi
				if [ ! -f "$SAM_DIR/aln${m}B.sam" -o $OVERWRITE_ALIGNMENTS -ge 1 ]; then
					echo "bwa mem -t $NCPU $REF_DIR/refB.fsa $READ_DIR/sim${m}.1.fq $READ_DIR/sim${m}.2.fq > $SAM_DIR/aln${m}B.sam"
					if [ $TEST == "0" ]; then
						#wait_on_cpu
						bwa mem -t $NCPU $REF_DIR/refB.fsa $READ_DIR/sim${m}.1.fq $READ_DIR/sim${m}.2.fq > $SAM_DIR/aln${m}B.sam
					fi
				fi
			done
		done
		#wait_on_finish
	done
done
