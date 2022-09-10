#!/bin/bash
# @author R. Kulkarni
# @file err2info.sh
#
# Purpose: part of the peant real data analysis pipeline
# convert stderr error output of capg_wgs to "info" format

TEST_ONLY=0
CAPG_HOME=.
INDIR=${CAPG_HOME}/data/err
OUTDIR=${CAPG_HOME}/data/info
CMD=${CAPG_HOME}/script_analysis/err2info.py

if [ ! -d $OUTDIR ]; then
	echo "mkdir --parents $OUTDIR"
	if [ $TEST_ONLY -ge 1 ]
	then
		mkdir --parents $OUTDIR
	fi
fi

for i in $INDIR/*.err;
do
	file=`basename $i .err`
	echo $file
	if [ "$i" -nt "$OUTDIR/${file}_info.txt" ]; then
		#echo "python $CMD -i $i -o $OUTDIR/${file}_info.txt"
		if [ $TEST_ONLY -ge 1 ]; then
			echo "python $CMD -i $i -o $OUTDIR/${file}_info.txt -e"
		else
			python $CMD -i $i -o $OUTDIR/${file}_info.txt -e	# with filters, coverage defaults to >=8
			#python $CMD -i $i -o $OUTDIR/${file}_info.txt -c 0	# without filters
		fi
		#python err2info_file.py -i $i -o $OUTDIR/${file}_info.txt
	fi
done

