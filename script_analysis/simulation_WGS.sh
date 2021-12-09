# This script generates simulation for WGS
#INPUTS:Reference sequence of one subgenome and error profile files
#OUTPUTS: Simulated reads, alignments and sam files 

#!/bin/bash

REFERENCE=/Path/ref.fa
OUTPUT_SIMU=/Path to output 
ERR_FILE1=/Path/ERRR1.txt
ERR_FILE2=/Path/ERRR2.txt

# Simulating under different homeologous rates
for p in 0.005 0.007 0.010
do
	PARENT_FOLDER="$OUTPUT_SIMU/homr${p}"
	mkdir "$PARENT_FOLDER"
	CMD="/home/peanut/peanut2/WGS/simulation_V6/ksd-roshan_wgs/run_simu -e $REFERENCE -j ${p} -g 0.005 -a 100 -b 100 -s 1 -o $PARENT_FOLDER/indiv -f $PARENT_FOLDER/ref -r $PARENT_FOLDER/ref.sam -n 50" 
	$CMD 2> $PARENT_FOLDER/truth_${p}.txt &
	bwa index $PARENT_FOLDER/refA.fsa
	bwa index $PARENT_FOLDER/refB.fsa
	mkdir "$PARENT_FOLDER/pair"
        # Simulating under differnt sequence coverage
	for q in 5  15 20
	do
		DATA_FOLDER="$PARENT_FOLDER/pair/cov${q}"
		mkdir "$DATA_FOLDER"
		for m in {0..49}
		do
			art_illumina -1 $ERR_FILE1 -2 $ERR_FILE2 -l 150 -i $PARENT_FOLDER/indiv${m}.fsa -o $DATA_FOLDER/sim${m} -f ${q} -p -m 300 -s 10
			bwa mem $PARENT_FOLDER/refA.fsa $DATA_FOLDER/sim${m}1.fq $DATA_FOLDER/sim${m}2.fq > $DATA_FOLDER/aln${m}A.sam
			bwa mem $PARENT_FOLDER/refB.fsa $DATA_FOLDER/sim${m}1.fq $DATA_FOLDER/sim${m}2.fq > $DATA_FOLDER/aln${m}B.sam
		done
	done
done
		

 
