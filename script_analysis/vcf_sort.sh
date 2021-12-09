IN_1=/path/indv_vcf_A
IN_2=/path/indv_vcf_B
OUT_1=/path/indv_vcf_A_sorted
OUT_2=/path/indv_vcf_B_sorted  

for GENOTYPE in SRR4124062 SRR4124066 SRR4124068 SRR4124074 SRR4124078 SRR8361734 SRR8361735 SRR8361736 SRR8361737 SRR8361738 SRR8736998 SRR8737008 SRR8737061 SRR8737062;
do
        for TARGET in $(seq 1 1000)
        do
		java -jar picard.jar SortVcf I=$IN_1/target_$TARGET/${GENOTYPE}_${TARGET}-A.vcf SEQUENCE_DICTIONARY=/home/peanut/peanut2/WGS/WGS_genotyping/new/tet_A.dict  O=$OUT_1/target_$TARGET/${GENOTYPE}_${TARGET}-A.vcf
	        java -jar picard.jar SortVcf I=$IN_2/${GENOTYPE}_${TARGET}-B.vcf SEQUENCE_DICTIONARY=/home/peanut/peanut2/WGS/WGS_genotyping/new/tet_B.dict  O=$OUT_2/target_$TARGET/${GENOTYPE}_${TARGET}-B.vcf
	done
done	
