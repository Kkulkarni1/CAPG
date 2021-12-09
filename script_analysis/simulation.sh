#!/bin/bash
#
# Purpose: Simulate WGS data.
#
# Input:
# $REFERENCE				Arahy.chr01:1-100000 from Tifrunner assembly, peanut subgenome A reference	
# $ERR_FILE[12]				miseq error rates			
# $SIMULATOR				Yudi's simulator
#
# Output:
# $PARENT_FOLDER/homr*/ref[AB].fsa	subgenomic references: refA.fsa identical to $REFERENCE, refB.fsa simulated
# $PARENT_FOLDER/homr*/ref.sam		alignment of ref[AB].fsa
# $PARENT_FOLDER/homr*/truth_*.txt	captured output of simulator showing all SNPs; loci are 0-based indices
# $PARENT_FOLDER/homr*/indiv*.fsa	individual's realized four genomes
# $PARENT_FOLDER/homr*/pair/cov*/aln*.sam	bwa aligned simulated reads	

SIMULATOR="~/scripts/ksd-roshan_wgs/run_simu"	# simulator executable		
OUTPUT_SIMU=/home/peanut/peanut2/CAPG/WGS/simulation				# output directory for simulation
REFERENCE="$OUTPUT_SIMU/ref_simulation.fa"					# reference for simulation	
ERR_FILE1="$OUTPUT_SIMU/miseq250R1.txt"						# art_illumina error files	
ERR_FILE2="$OUTPUT_SIMU/miseq250R2.txt"

# simulation conditions
ALPHA=100		# SIMULATOR: allele proportions from Beta(ALPHA, BETA)
BETA=100
ALLELIC_RATE=0.005	# SIMULATOR: probability of allelic SNP at sites that are not homoeologous loci
NREP=50			# SIMULATOR: number of genotypes to simulate
READ_LEN=150		# ART: read length
FRAG_LEN_MEAN=300	# ART: fragment length mean
FRAG_LEN_SD=10		# ART: fragment length standard deviation

EXE=roshan
NCPU=7

function wait_on {
	NUM=`pgrep -c $EXE`
	while [ $NUM -ge $NCPU ]
	do
		echo "$NUM $EXE running..."
		sleep 2
		NUM=`pgrep -c $EXE`
	done
}

for p in 0.005 0.007 0.010	# homoeologous rates
do
	PARENT_FOLDER="$OUTPUT_SIMU/homr${p}"
	mkdir "$PARENT_FOLDER"
	CMD="$SIMULATOR -e $REFERENCE -j ${p} -g $ALLELIC_RATE -a $ALPHA -b $BETA -s 1 -o $PARENT_FOLDER/indiv -f $PARENT_FOLDER/ref -r $PARENT_FOLDER/ref.sam -n $NREP"
	wait_on
        # Simulating reference B genome, allotetraploid individuals and generating truth files
	$CMD 2> $PARENT_FOLDER/truth_${p}.txt &
	bwa index $PARENT_FOLDER/refA.fsa
	bwa index $PARENT_FOLDER/refB.fsa
	mkdir "$PARENT_FOLDER/pair"
	for q in 5 10 15 20	# coverage levels
	do
		DATA_FOLDER="$PARENT_FOLDER/pair/cov${q}"
		mkdir "$DATA_FOLDER"
		for m in {0..49} #number of individuals
		do
                        # Simulating reads
			art_illumina -1 $ERR_FILE1 -2 $ERR_FILE2 -l $READ_LEN -i $PARENT_FOLDER/indiv${m}.fsa -o $DATA_FOLDER/sim${m} -f ${q} -p -m $FRAG_LEN_MEAN -s $FRAG_LEN_SD
                        # Aligning reads to A and B subgenomes
			bwa mem $PARENT_FOLDER/refA.fsa $DATA_FOLDER/sim${m}1.fq $DATA_FOLDER/sim${m}2.fq > $DATA_FOLDER/aln${m}A.sam
			bwa mem $PARENT_FOLDER/refB.fsa $DATA_FOLDER/sim${m}1.fq $DATA_FOLDER/sim${m}2.fq > $DATA_FOLDER/aln${m}B.sam
		done
	done
done
