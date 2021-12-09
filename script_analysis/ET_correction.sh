IN_1=/path/merged_vcfA
IN_2=/path/merged_vcfB
OUT_1=/path/merged_vcfA_ET_filter
OUT_2=/path/merged_vcfB_ET_filter

for TARGET in `seq 0 1000`
do
	python ET_correction.py -i $IN_1/merged_${TARGET}.txt -o $OUT_1/merged_${TARGET}_A_ET.txt
	python ET_correction.py -i $IN_2/merged_${TARGET}.txt -o $OUT_2/merged_${TARGET}_B_ET.txt
done
