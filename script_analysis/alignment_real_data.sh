# This script aligns FASTQ reads using bwa
# Input: FASTQ files
# Output: bamfiles

INDIR=/path/fastqfiles
OUTDIR=/path/samfiles

for i in INDIR;
do
	file=`basename $i _1.fastq`
	bwa-mem2 -t32 /path/refA.fsa /path/${file}_1.fastq /path/${file}_2.fastq | samtools sort -o /path/${file}_A.bam
	bwa-mem2 -t32 /path/refB.fsa /path/${file}_1.fastq /path/${file}_2.fastq | samtools sort -o /path/${file}_B.bam
done  
