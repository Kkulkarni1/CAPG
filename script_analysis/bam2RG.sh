# This script adds read groups to bamfiles 
# Input: bamfiles without read group
# Output: bamfiles with read group

IN_1=/path/BAMA_target
IN_2=/path/BAMB_target
OUT_1=/path/BAMA_target_RG
OUT_2=/path/BAMB_target_RG

for GENOTYPE in SRR4124062 SRR4124066 SRR4124068 SRR4124074 SRR4124078 SRR4124079 SRR8361734 SRR8361735 SRR8361736 SRR8361737 SRR8361738 SRR8736998 SRR8737001 SRR8737008 SRR8737061 SRR8737062
do
	java -jar /opt/picard/picard.jar AddOrReplaceReadGroups I=$IN_1/${GENOTYPE}_A.target.bam O=$OUT_1/${GENOTYPE}_A.target.bam RGID=${GENOTYPE} RGPL=illumina RGPU=uni1 RGSM=${GENOTYPE} RGLB=lib1

	java -jar /opt/picard/picard.jar AddOrReplaceReadGroups I=$IN_2/${GENOTYPE}_B.target.bam O=$OUT_2/${GENOTYPE}_B.target.bam RGID=${GENOTYPE} RGPL=illumina RGPU=uni1 RGSM=${GENOTYPE} RGLB=lib1

done

	
