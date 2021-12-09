#!/bin/bash
# This script performs genotyping and SNP calling on simulation data using GATK
# INPUT: conacatenated reference file, A and B genome reference files and simukated fastq reads
# OUTPUT: vcf files
 
DATA_PATH=/path

for p in 0.005 0.007 0.010
for p in 0.005
do
	PARENT_FOLDER="$DATA_PATH/homr${p}"
	cat $PARENT_FOLDER/refA.fa $PARENT_FOLDER/refB.fa > $PARENT_FOLDER/ref_all.fa
	bwa index $PARENT_FOLDER/ref_all.fa

	# preparing index files for alignment
	bwa index $PARENT_FOLDER/refA.fa
	bwa index $PARENT_FOLDER/refB.fa
	samtools faidx $PARENT_FOLDER/refA.fa
	samtools faidx $PARENT_FOLDER/refB.fa
	
	# prpeparing dictionary file needed by GATK
	java -jar /usr/local/jar/picard.jar CreateSequenceDictionary REFERENCE=$PARENT_FOLDER/refA.fa OUTPUT=$PARENT_FOLDER/refA.dict
	java -jar /usr/local/jar/picard.jar CreateSequenceDictionary REFERENCE=$PARENT_FOLDER/refB.fa OUTPUT=$PARENT_FOLDER/refB.dict

	for q in 5 10 20
	do
		DATA_FOLDER="$PARENT_FOLDER/pair/cov${q}"

		for m in {0..49}
		do
			# split the reads to A and B
			bwa mem $PARENT_FOLDER/ref_all.fa $DATA_FOLDER/sim_files/sim${m}1.fq $DATA_FOLDER/sim_files/sim${m}2.fq > $DATA_FOLDER/GATK/sim${m}test.sam
			samtools view -S -b $DATA_FOLDER/GATK/sim${m}test.sam > $DATA_FOLDER/GATK/sim${m}test.bam
			samtools sort $DATA_FOLDER/GATK/sim${m}test.bam -o $DATA_FOLDER/GATK/sim${m}index.test.bam
			samtools index $DATA_FOLDER/GATK/sim${m}index.test.bam
			samtools view -b $DATA_FOLDER/GATK/sim${m}index.test.bam Genome_B > $DATA_FOLDER/GATK/sim${m}in_chr1.bam
			samtools view -b $DATA_FOLDER/GATK/sim${m}index.test.bam Genome_A > $DATA_FOLDER/GATK/sim${m}in_chr0.bam
			bamToFastq -i $DATA_FOLDER/GATK/sim${m}in_chr1.bam -fq $DATA_FOLDER/GATK/sim${m}B.fastq
			bamToFastq -i $DATA_FOLDER/GATK/sim${m}in_chr0.bam -fq $DATA_FOLDER/GATK/sim${m}A.fastq

			bwa mem $PARENT_FOLDER/refA.fa $DATA_FOLDER/GATK/sim${m}A.fastq > $DATA_FOLDER/GATK/sim${m}A.sam
			bwa mem $PARENT_FOLDER/refB.fa $DATA_FOLDER/GATK/sim${m}B.fastq > $DATA_FOLDER/GATK/sim${m}B.sam

			# sam to bam conversion
			samtools view -bS $DATA_FOLDER/GATK/sim${m}A.sam > $DATA_FOLDER/GATK/sim${m}A.bam
			samtools view -bS $DATA_FOLDER/GATK/sim${m}B.sam > $DATA_FOLDER/GATK/sim${m}B.bam

			# Adding @RG to header
			java -jar /usr/local/jar/picard.jar AddOrReplaceReadGroups I=$DATA_FOLDER/GATK/other_files/sim${m}A.bam O=$DATA_FOLDER/GATK/other_files/sim${m}A-RG.bam RGID=aln$m RGPL=illumina RGPU=uni1 RGSM=aln$m RGLB=lib1
			java -jar /usr/local/jar/picard.jar AddOrReplaceReadGroups I=$DATA_FOLDER/GATK/other_files/sim${m}B.bam O=$DATA_FOLDER/GATK/other_files/sim${m}B-RG.bam RGID=aln$m RGPL=illumina RGPU=uni1 RGSM=aln$m RGLB=lib1

			# Sort bam files
			samtools sort $DATA_FOLDER/GATK/other_files/sim${m}A-RG.bam -o $DATA_FOLDER/GATK/other_files/sim${m}A.sorted.bam
			samtools sort $DATA_FOLDER/GATK/other_files/sim${m}B-RG.bam -o $DATA_FOLDER/GATK/other_files/sim${m}B.sorted.bam

			# index sorted bam files
			samtools index $DATA_FOLDER/GATK/other_files/sim${m}A.sorted.bam
			samtools index $DATA_FOLDER/GATK/other_files/sim${m}B.sorted.bam

			# Call SNPs using GATK
			/usr/local/bin/gatk HaplotypeCaller -R $PARENT_FOLDER/refA.fa -I $DATA_FOLDER/GATK/other_files/sim${m}A.sorted.bam -O $DATA_FOLDER/GATK/sim${m}A.vcf -ERC BP_RESOLUTION
			/usr/local/bin/gatk HaplotypeCaller -R $PARENT_FOLDER/refB.fa -I $DATA_FOLDER/GATK/other_files/sim${m}B.sorted.bam -O $DATA_FOLDER/GATK/sim${m}B.vcf -ERC BP_RESOLUTION
			done
	done
done


	
