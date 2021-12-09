## For Simulation ##

INDIR_1=/home/peanut/peanut2/CAPG/WGS/simulation/homr0.010/pair/cov20/err_files_new 
#OUTDIR_1=/home/peanut/peanut2/CAPG/WGS/simulation/homr0.010/pair/cov20/info_files_new
OUTDIR_1=/home/peanut/peanut2/CAPG/WGS/simulation/homr0.010/pair/cov20/info_files_new_w_filter

for i in $INDIR_1/*.err;
do
	file=`basename $i .err`
	python err2info_file.py -i $i -o $OUTDIR_1/${file}_info.txt
done
