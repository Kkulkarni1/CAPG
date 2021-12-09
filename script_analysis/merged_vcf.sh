IN_1=/path/vcfA_filt_cov
IN_2=/path/vcfB_filt_cov
OUT_1=/path/merged_vcfA
OUT_2=/path/merged_vcfB

for GENOTYPE in SRR4124062 SRR4124066 SRR4124068 SRR4124074 SRR4124078 SRR8361734 SRR8361735 SRR8361736 SRR8361737 SRR8361738 SRR8736998 SRR8737008 SRR8737061 SRR8737062;
do
	for TARGET in $(seq 1 1000)
	do
		/home/peanut/peanut2/CAPG/WGS/real_data/peanut_V2/htslib/bgzip -f $IN_1/target_${TARGET}/${GENOTYPE}_${TARGET}-A.vcf
		/home/peanut/peanut2/CAPG/WGS/real_data/peanut_V2/htslib/tabix -f $IN_1/target_${TARGET}/${GENOTYPE}_${TARGET}-A.vcf.gz
		/home/peanut/peanut2/CAPG/WGS/real_data/peanut_V2/htslib/bgzip -f $IN_2/target_${TARGET}/${GENOTYPE}_${TARGET}-B.vcf
		/home/peanut/peanut2/CAPG/WGS/real_data/peanut_V2/htslib/tabix -f $IN_2/target_${TARGET}/${GENOTYPE}_${TARGET}-B.vcf.gz
	done
done

for TARGET in $(seq 1 1000)
do
        bcftools merge $IN_1/target_${TARGET}/*.gz > $OUT_1/merged_${TARGET}.txt
	bcftools merge $IN_2/target_${TARGET}/*.gz > $OUT_2/merged_${TARGET}.txt
done
