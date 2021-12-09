# This script downloads FASTQ files based on SRA ID
# Input: SRA.list file that contains SRA IDs

cat SRA.list | while read line;
do
	fasterq-dump --split-files ${line}
done
