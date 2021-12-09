IN_DIR=/Path/fasta_files
OUT_DIR=/Path/OUTDIR

for i in $IN_DIR/*.fsa;
do
	file=`basename $i .fsa`;
	Rscript truth.R -i $i  --out $OUT_DIR/${file}_truth_het_homeo.txt
done
