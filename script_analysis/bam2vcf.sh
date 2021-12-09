# This script calls SNPs using GATK's HaplotypeCaller
# Input: sorted bam files
# Output: vcf files

IN_1=/path/BAMA_target_RG
IN_2=/path/BAMB_target_RG
OUT_1=/path/vcf_A
OUT_2=/path/vcf_B

for GENOTYPE in SRR4124062 SRR4124066 SRR4124068 SRR4124074 SRR4124078 SRR4124079 SRR8361734 SRR8361735 SRR8361736 SRR8361737 SRR8361738 SRR8736998 SRR8737001 SRR8737008 SRR8737061 SRR8737062
do

        /path/gatk HaplotypeCaller -R target_ref_gap_A.fasta  -I $IN_1/${GENOTYPE}_A.target.bam -O $OUT_1/${GENOTYPE}_A.vcf -ERC BP_RESOLUTION
        /path/gatk HaplotypeCaller -R target_ref_gap_B.fasta  -I $IN_2/${GENOTYPE}_B.target.bam -O $OUT_2/${GENOTYPE}_B.vcf -ERC BP_RESOLUTION
done    
