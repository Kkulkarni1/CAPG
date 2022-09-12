#!/usr/bin/bash
# @file err2info_simulation.sh
#
# Purpose: convert CAPG captured stderr output to "info" files

TEST=0					# test without doing anything
OVERWRITE=1				# overwrite "info" files
CAPG_HOME=.				# top of CAPG repository
SIM_DIR=$CAPG_HOME/data/simulation	# IN/OUTPUT: simulation data
EXE_DIR=script_analysis			# INPUT: where the python script is
ERR_DIR=err				# INPUT: captured CAPG stderr output files
INFO_DIR=info				# OUTPUT: where to put "info" files
HR_RATES="0.005 0.007"			# SIMULATOR: homoeologous SNP rate (manuscript: "0.005 0.007 0.100")
COV_RATES="5 10"			# SIMULATOR: coverage rates per chromosome (manuscript: "5 10 20" for hr = 0.007)
MM_RATES="0.000 0.001 0.010"		# SIMULATOR: mismatch rate (manuscript w/ mistake: "0.00 0.00 0.01" labeled as "0.00 0.01 0.10")
NREP=1					# SIMULATOR: number of genotypes to simulate (50)

INDIR=/home/peanut/peanut2/CAPG/WGS/simulation/homr0.007mm0.010/pair/cov5/err_files
OUTDIR=/home/peanut/peanut2/CAPG/WGS/simulation/homr0.007mm0.010/pair/cov5/info_files


for hr in $HR_RATES
do
	for cov in $COV_RATES
	do
		for mm in $MM_RATES
		do
			DIR="$SIM_DIR/homr${hr}mm${mm}/cov${cov}"

			if [ ! -d "$DIR/$INFO_DIR" ]; then
				echo "mkdir --parents $DIR/$INFO_DIR"
				if [ $TEST == "0" ]; then
					mkdir --parents $DIR/$INFO_DIR
				fi
			fi

			for i in $DIR/$ERR_DIR/*.err
			do
				file=`basename $i .err`
				if [ ! -f "$DIR/$INFO_DIR/${file}_info.txt" -o $OVERWRITE -ge 1 ]; then
					echo "python $EXE_DIR/err2info_simulation.py --in_file $i --out_file $DIR/$INFO_DIR/${file}_info.txt"
					if [ $TEST = "0" ]; then
						python $EXE_DIR/err2info_simulation.py --in_file $i --out_file $DIR/$INFO_DIR/${file}_info.txt
					fi
				fi
			done
		done
	done
done
