# This script subsets bam files based on coordinates in the bed file
# Input: bam files, target bed files

INDIR=/path/bamfiles/
OUTDIR=/path/bam_subset/

for i in bamfiles
do
	file=`basename $i _A.bam`;
	samtools view -b -L target_genes_A.bed $INDIR/${file}_A.bam > $OUTDIR/${file}_A_subset.bam
	samtools view -b -L target_genes_B.bed $INDIR/${file}_B.bam > $OUTDIR/${file}_B_subset.bam

done 
