#!/usr/bin/bash
# @file prep_gatk_simulation.sh
# @author R. Kulkarni
#
# Requires: gatk, samtools, bwa
#
# This script performs genotyping and SNP calling on simulation data using GATK
# INPUT: conacatenated reference file, A and B genome reference files and simukated fastq reads
# OUTPUT: vcf files
 
TEST=0				# test without doing anything
OVERWRITE=1			# overwrite vcf files
CAPG_HOME=.			# top of CAPG repository
GATK=/opt/gatk-4.2.6.1/gatk	# gatk executable
SIM_DIR=data/simulation		# IN/OUTPUT: simulation data
VCF_DIR=vcf			# OUTPUT: vcf directory
REF_SAM=ref.sam			# INPUT: references, including location specification
REF_FSA=ref.fsa			# INPUT: name of reference FASTA file
REFA_FSA=refA.fsa		# INPUT: name of subgenome A reference FASTA file
REFB_FSA=refB.fsa		# INPUT: name of subgenome B reference FASTA file
HR_RATES="0.005 0.007"		# SIMULATOR: homoeologous SNP rate (manuscript: "0.005 0.007 0.100")
COV_RATES="5 10"		# SIMULATOR: coverage rates per chromosome (manuscript: "5 10 20" for hr = 0.007)
MM_RATES="0.000 0.001 0.010"	# SIMULATOR: mismatch rate (manuscript w/ mistake: "0.00 0.00 0.01" labeled as "0.00 0.01 0.10")
NREP=1				# SIMULATOR: number of genotypes to simulate (50)
MM_BASE="0.000"			# identify where the FASTQ read files are


for hr in $HR_RATES; do
for mm in $MM_RATES; do
	REF_DIR="$CAPG_HOME/$SIM_DIR/homr${hr}mm${mm}"

	REF_NAME_A=`head -n 1 $REF_DIR/$REF_SAM | cut -f 2 | sed 's/SN://'`
	REF_NAME_B=`tail -n +2 $REF_DIR/$REF_SAM | cut -f 1`

	if [ ! -f "$REF_DIR/$REF_FSA.bwt" -o $OVERWRITE -ge 1 ]; then
		echo "bwa index $REF_DIR/$REF_FSA"
		if [ $TEST == "0" ]; then
			bwa index $REF_DIR/$REF_FSA
		fi
	fi

	# preparing index files for alignment
	if [ ! -f "$REF_DIR/$REFA_FSA.bwt" -o $OVERWRITE -ge 1 ]; then
		echo "bwa index $REF_DIR/$REFA_FSA"
		if [ $TEST == "0" ]; then
			bwa index $REF_DIR/$REFA_FSA
		fi
	fi
	if [ ! -f "$REF_DIR/$REFB_FSA.bwt" -o $OVERWRITE -ge 1 ]; then
		echo "bwa index $REF_DIR/$REFB_FSA"
		if [ $TEST == "0" ]; then
			bwa index $REF_DIR/$REFB_FSA
		fi
	fi
	
	# prpeparing dictionary file needed by GATK
	if [ ! -f "$REF_DIR/refA.dict" -o $OVERWRITE -ge 1 ]; then
		# stupid gatk restriction
		if [ ! -f "$REF_DIR/refA.fasta" -o $OVERWRITE -ge 1 ]; then
			echo "ln -s -f $REFA_FSA $REF_DIR/refA.fasta"
			if [ $TEST == "0" ]; then
				ln -s -f $REFA_FSA $REF_DIR/refA.fasta
			fi
		fi
		if [ -f "$REF_DIR/refA.dict" ]; then
			echo "rm $REF_DIR/refA.dict"
			rm $REF_DIR/refA.dict
		fi
		echo "$GATK CreateSequenceDictionary -REFERENCE $REF_DIR/refA.fasta -OUTPUT $REF_DIR/refA.dict"
		if [ $TEST == "0" ]; then
			$GATK CreateSequenceDictionary -REFERENCE $REF_DIR/refA.fasta -OUTPUT $REF_DIR/refA.dict
		fi
	fi
	if [ ! -f "$REF_DIR/refB.dict" -o $OVERWRITE -ge 1 ]; then
		# stupid gatk restriction
		if [ ! -f "$REF_DIR/refB.fasta" -o $OVERWRITE -ge 1 ]; then
			echo "ln -s -f $REFB_FSA $REF_DIR/refB.fasta"
			if [ $TEST == "0" ]; then
				ln -s -f $REFB_FSA $REF_DIR/refB.fasta
			fi
		fi
		if [ -f "$REF_DIR/refB.dict" ]; then
			echo "rm $REF_DIR/refB.dict"
			rm $REF_DIR/refB.dict
		fi
		echo "$GATK CreateSequenceDictionary -REFERENCE $REF_DIR/refB.fasta -OUTPUT $REF_DIR/refB.dict"
		if [ $TEST == "0" ]; then
			$GATK CreateSequenceDictionary -REFERENCE $REF_DIR/refB.fasta -OUTPUT $REF_DIR/refB.dict
		fi
	fi
	if [ ! -f "$REF_DIR/refA.fasta.fai" -o $OVERWRITE -ge 1 ]; then
		echo "samtools faidx $REF_DIR/refA.fasta"
		if [ $TEST == "0" ]; then
			samtools faidx $REF_DIR/refA.fasta
		fi
	fi
	if [ ! -f "$REF_DIR/refB.fasta.fai" -o $OVERWRITE -ge 1 ]; then
		echo "samtools faidx $REF_DIR/refB.fasta"
		if [ $TEST == "0" ]; then
			samtools faidx $REF_DIR/refB.fasta
		fi
	fi

	for cvg in $COV_RATES; do
		READ_DIR="$CAPG_HOME/$SIM_DIR/homr${hr}mm$MM_BASE/cov$cvg"
		OUT_DIR="$REF_DIR/cov$cvg"

		if [ ! -d "$REF_DIR/GATK/bam" ]; then
			echo "mkdir --parents $OUT_DIR/GATK/bam"
			if [ $TEST == "0" ]; then
				mkdir --parents $OUT_DIR/GATK/bam
			fi
		fi
		if [ ! -d "$REF_DIR/GATK/$VCF_DIR" ]; then
			echo "mkdir --parents $OUT_DIR/GATK/$VCF_DIR"
			if [ $TEST == "0" ]; then
				mkdir --parents $OUT_DIR/GATK/$VCF_DIR
			fi
		fi


		for m in `seq 1 $NREP`; do
			m=$((m-1))

			# partition the reads to A and B
			# this is what we did, even though it is probably overkill
			if [ ! -f "$OUT_DIR/GATK/bam/sim${m}.sorted.bam" -o $OVERWRITE -ge 1 ]; then
				echo "bwa mem $REF_DIR/$REF_FSA $READ_DIR/fastq/sim${m}.1.fq $READ_DIR/fastq/sim${m}.2.fq > $OUT_DIR/GATK/bam/sim${m}.all.sam"
				echo "samtools view -S -b $OUT_DIR/GATK/bam/sim${m}.sam > $OUT_DIR/GATK/bam/sim${m}.all.bam"
				echo "samtools sort $OUT_DIR/GATK/bam/sim${m}.bam -o $OUT_DIR/GATK/bam/sim${m}.sorted.all.bam"
				echo "samtools index $OUT_DIR/GATK/bam/sim${m}.sorted.all.bam"
				if [ $TEST == "0" ]; then
					bwa mem $REF_DIR/$REF_FSA $READ_DIR/fastq/sim${m}.1.fq $READ_DIR/fastq/sim${m}.2.fq > $OUT_DIR/GATK/bam/sim${m}.all.sam
					samtools view -S -b $OUT_DIR/GATK/bam/sim${m}.sam > $OUT_DIR/GATK/bam/sim${m}.all.bam
					samtools sort $OUT_DIR/GATK/bam/sim${m}.bam -o $OUT_DIR/GATK/bam/sim${m}.sorted.all.bam
					samtools index $OUT_DIR/GATK/bam/sim${m}.sorted.all.bam
				fi
			fi

			if [ ! -f "$OUT_DIR/GATK/bam/sim${m}A.bam -o ! -f "$OUT_DIR/GATK/bam/sim${m}B.bam -o $OVERWRITE -ge 1 ]; then
				echo "samtools view -b $OUT_DIR/GATK/bam/sim${m}.sorted.all.bam $REF_NAME_A > $OUT_DIR/GATK/bam/sim${m}A.split.bam"
				echo "samtools view -b $OUT_DIR/GATK/bam/sim${m}.sorted.all.bam $REF_NAME_B > $OUT_DIR/GATK/bam/sim${m}B.split.bam"
				echo "bamToFastq -i $OUT_DIR/GATK/bam/sim${m}A.split.bam -fq $OUT_DIR/GATK/bam/sim${m}A.fq"
				echo "bamToFastq -i $OUT_DIR/GATK/bam/sim${m}B.split.bam -fq $OUT_DIR/GATK/bam/sim${m}B.fq"
				echo "bwa mem $REF_DIR/$REFA_FSA $OUT_DIR/GATK/bam/sim${m}A.fq > $OUT_DIR/GATK/bam/sim${m}A.sam"	# overwrite
				echo "bwa mem $REF_DIR/$REFB_FSA $OUT_DIR/GATK/bam/sim${m}B.fq > $OUT_DIR/GATK/bam/sim${m}B.sam"	# overwrite
				echo "samtools view -bS $OUT_DIR/GATK/bam/sim${m}A.sam > $OUT_DIR/GATK/bam/sim${m}A.bam"	# overwrite
				echo "samtools view -bS $OUT_DIR/GATK/bam/sim${m}B.sam > $OUT_DIR/GATK/bam/sim${m}B.bam"	# overwrite
				if [ $TEST == "0" ]; then
					samtools view -b $OUT_DIR/GATK/bam/sim${m}.sorted.all.bam $REF_NAME_A > $OUT_DIR/GATK/bam/sim${m}A.split.bam
					samtools view -b $OUT_DIR/GATK/bam/sim${m}.sorted.all.bam $REF_NAME_A > $OUT_DIR/GATK/bam/sim${m}B.split.bam
					bamToFastq -i $OUT_DIR/GATK/bam/sim${m}A.split.bam -fq $OUT_DIR/GATK/bam/sim${m}A.fq
					bamToFastq -i $OUT_DIR/GATK/bam/sim${m}B.split.bam -fq $OUT_DIR/GATK/bam/sim${m}B.fq
					bwa mem $REF_DIR/$REFA_FSA $OUT_DIR/GATK/bam/sim${m}A.fq > $OUT_DIR/GATK/bam/sim${m}A.sam	# overwrite
					bwa mem $REF_DIR/$REFB_FSA $OUT_DIR/GATK/bam/sim${m}B.fq > $OUT_DIR/GATK/bam/sim${m}B.sam	# overwrite
					samtools view -bS $OUT_DIR/GATK/bam/sim${m}A.sam > $OUT_DIR/GATK/bam/sim${m}A.bam	# overwrite
					samtools view -bS $OUT_DIR/GATK/bam/sim${m}B.sam > $OUT_DIR/GATK/bam/sim${m}B.bam	# overwrite
				fi
			fi

			# Adding @RG to header, sort bam files, and index sorted bam files
			if [ ! -f "$OUT_DIR/GATK/bam/sim${m}A.sorted.bam" -o $OVERWRITE -ge 1 ]; then
				echo "$GATK AddOrReplaceReadGroups -I $OUT_DIR/GATK/bam/sim${m}A.bam -O $OUT_DIR/GATK/bam/sim${m}A-RG.bam -RGID aln$m -RGPL illumina -RGPU uni1 -RGSM aln$m -RGLB lib1"
				echo "samtools sort $OUT_DIR/GATK/bam/sim${m}A-RG.bam -o $OUT_DIR/GATK/bam/sim${m}A.sorted.bam"
				echo "samtools index $OUT_DIR/GATK/bam/sim${m}A.sorted.bam"
				if [ $TEST == "0" ]; then
					$GATK AddOrReplaceReadGroups -I $OUT_DIR/GATK/bam/sim${m}A.bam -O $OUT_DIR/GATK/bam/sim${m}A-RG.bam -RGID aln$m -RGPL illumina -RGPU uni1 -RGSM aln$m -RGLB lib1
					samtools sort $OUT_DIR/GATK/bam/sim${m}A-RG.bam -o $OUT_DIR/GATK/bam/sim${m}A.sorted.bam
					samtools index $OUT_DIR/GATK/bam/sim${m}A.sorted.bam
				fi
			fi
			if [ ! -f "$OUT_DIR/GATK/bam/sim${m}B.sorted.bam" -o $OVERWRITE -ge 1 ]; then
				echo "$GATK AddOrReplaceReadGroups -I $OUT_DIR/GATK/bam/sim${m}B.bam -O $OUT_DIR/GATK/bam/sim${m}B-RG.bam -RGID aln$m -RGPL illumina -RGPU uni1 -RGSM aln$m -RGLB lib1"
				echo "samtools sort $OUT_DIR/GATK/bam/sim${m}B-RG.bam -o $OUT_DIR/GATK/bam/sim${m}B.sorted.bam"
				echo "samtools index $OUT_DIR/GATK/bam/sim${m}B.sorted.bam"
				if [ $TEST == "0" ]; then
					$GATK AddOrReplaceReadGroups -I $OUT_DIR/GATK/bam/sim${m}B.bam -O $OUT_DIR/GATK/bam/sim${m}B-RG.bam -RGID aln$m -RGPL illumina -RGPU uni1 -RGSM aln$m -RGLB lib1
					samtools sort $OUT_DIR/GATK/bam/sim${m}B-RG.bam -o $OUT_DIR/GATK/bam/sim${m}B.sorted.bam
					samtools index $OUT_DIR/GATK/bam/sim${m}B.sorted.bam
				fi
			fi


			# Call SNPs using GATK
			if [ ! -f "$OUT_DIR/GATK/$VCF_DIR/sim${m}A.vcf" -o $OVERWRITE -ge 1 ]; then
				echo "$GATK HaplotypeCaller -R $REF_DIR/refA.fasta -I $OUT_DIR/GATK/bam/sim${m}A.sorted.bam -O $OUT_DIR/GATK/$VCF_DIR/sim${m}A.vcf -ERC BP_RESOLUTION"
				if [ $TEST == "0" ]; then
					$GATK HaplotypeCaller -R $REF_DIR/refA.fasta -I $OUT_DIR/GATK/bam/sim${m}A.sorted.bam -O $OUT_DIR/GATK/$VCF_DIR/sim${m}A.vcf -ERC BP_RESOLUTION
				fi
			fi
			if [ ! -f "$OUT_DIR/GATK/$VCF_DIR/sim${m}B.vcf" -o $OVERWRITE -ge 1 ]; then
				echo "$GATK HaplotypeCaller -R $REF_DIR/refB.fasta -I $OUT_DIR/GATK/bam/sim${m}B.sorted.bam -O $OUT_DIR/GATK/$VCF_DIR/sim${m}B.vcf -ERC BP_RESOLUTION"
				if [ $TEST == "0" ]; then
					$GATK HaplotypeCaller -R $REF_DIR/refB.fasta -I $OUT_DIR/GATK/bam/sim${m}B.sorted.bam -O $OUT_DIR/GATK/$VCF_DIR/sim${m}B.vcf -ERC BP_RESOLUTION
				fi
			fi
		done
	done
done
done
