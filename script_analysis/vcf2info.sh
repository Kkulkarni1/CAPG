#!/bin/bash
#
# Purpose:  Extract information for genotype calling from GATK vcf files.
#

TEST_ONLY=0
OVERWRITE=1
DIR=/path
EXE=vcf2infoksd.py
OSDIR=/path
EXTENSION=_GATK_info.txt

for hr in '0.005' '0.007' '0.010'; do
	for cvg in 5 10 20; do
		if [[ ! -d "$DIR/homr${hr}/pair/cov${cvg}" ]]; then
			continue
		fi
		for i in `seq 0 49`; do
			if [[ $OVERWRITE -eq 0 ]]
			then
				NLINES=`wc --lines $DIR/$OSDIR/homr${hr}/cov${cvg}/info_files/aln${i}${EXTENSION} | awk '{print $1}'`
				if [[ -f "$DIR/$OSDIR/homr${hr}/cov${cvg}/info_files/aln${i}${EXTENSION}" && $NLINES -eq 100001 ]]
				then
					continue
				fi
			fi
			echo "$DIR/homr${hr}/pair/cov${cvg}/GATK/vcf_A/sim${i}A.vcf"
			if [[ -f "$DIR/homr${hr}/pair/cov${cvg}/GATK/vcf_A/sim${i}A.vcf" && -f "$DIR/homr${hr}/pair/cov${cvg}/GATK/vcf_B/sim${i}B.vcf" ]]; then
				echo "mkdir --parents $DIR/$OSDIR/homr${hr}/cov${cvg}/info_files"
				if [ $TEST_ONLY == 0 ]
				then
					mkdir --parents $DIR/$OSDIR/homr${hr}/cov${cvg}/info_files
				fi
				echo "$EXE -i1 $DIR/homr${hr}/pair/cov${cvg}/GATK/vcf_A/sim${i}A.vcf -i2 $DIR/homr${hr}/pair/cov${cvg}/GATK/vcf_B/sim${i}B.vcf -o $DIR/$OSDIR/homr${hr}/cov${cvg}/info_files/aln${i}${EXTENSION}"
				if [ $TEST_ONLY == 0 ]
				then
					$EXE -i1 $DIR/homr${hr}/pair/cov${cvg}/GATK/vcf_A/sim${i}A.vcf -i2 $DIR/homr${hr}/pair/cov${cvg}/GATK/vcf_B/sim${i}B.vcf -o $DIR/$OSDIR/homr${hr}/cov${cvg}/info_files/aln${i}${EXTENSION}
				fi
			fi
		done
	done
done
