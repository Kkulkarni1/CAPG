# This scripts runs CAPG to perform genotyping 
# INPUTS: refernce files, SAM files, Alignment sam file
# OUPUTS: Info files with posterior probability, vcf files

REFA=/Path/refA.fa
REFB=/Path/refB.fa
SAMA=/Path/SAMA
SAMB=/Path/SAMB
EXTRACTED=/Path/extracted
VCFA=/Path/vcf_A
VCFB=/Path/vcf_B
ERR_FILE=/Path/err_files

for GENOTYPE in aln0 aln1 aln2 aln3 aln4 aln5 aln6 aln7 aln8 aln9 aln10 aln11 aln12 aln13 aln14 aln15 aln16 aln17 aln18 aln19 aln20 aln21 aln22 aln23 aln24 aln25 aln26 aln27 aln28 aln29 aln30 aln31 aln32 aln33 aln34 aln35 aln36 aln37 aln38 aln39 aln40 aln41 aln42 aln43 aln44 aln45 aln46 aln47 aln48 aln49;
do
	#file=`basename $i A.sam`;
	CMD="/Path/roshan --geno /Path/ref.sam --fsa_files $REFA $REFB --sam_files $SAMA/${GENOTYPE}A.sam $SAMB/${GENOTYPE}B.sam --ref_names Genome_A:0-100000 Genome_B:0-100000 -j $EXTRACTED/${GENOTYPE}_extracted.fsa --vcf_files $VCFA/${GENOTYPE}-A.vcf $VCFB/${GENOTYPE}-B.vcf --name ${GENOTYPE}"
	wait_on
	$CMD 2> $ERR_FILE/${GENOTYPE}.err &
done
